#include "regularity_check.h"
#include "input.h"
#include "poset.h"

#include <CGAL/number_utils.h> 

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>



// for each node in P1 we reconstruct CGAL triangulation with replay function
//export it (points + triangles) into Macaulay2 script
// run isRegularTriangulation and parse the true/false result 
// we print a summary table

// fix: do not run Macaulay2 for every triangulation but maybe pass all triangulations in one script 
// no urgent fix, just improvement 

// we use WL (windows subsystem for linux) to run M2 
// so we have to convert windows paths into WSL paths and run M2 via WSL command line

namespace regcheck {

namespace {

// we create a temporary M2 script for each triangulation
// they all need different names so we do not overwrite them 
// they are named regcheck_0.m2, regcheck_1.m2, ... and are deleted after running
// ---------- temp filename ----------
static std::string tmp_filename(const std::string& prefix, const std::string& ext) {
    static int counter = 0;
    std::ostringstream oss;
    oss << prefix << "_" << counter++ << ext;
    return oss.str();
}


// since we have to use WSL to run M2, we need to convert windows path into WSL path
// so convert windows path into linux path for WSL: C:\path\to\file -> /mnt/c/path/to/file
// ---------- Windows path -> WSL (/mnt/c/...) ----------
static std::string win_to_wsl_path(std::string win) {
    // convert backslashes
    std::replace(win.begin(), win.end(), '\\', '/');

    // if looks like "C:/..."
    if (win.size() >= 3 && std::isalpha((unsigned char)win[0]) && win[1] == ':' && win[2] == '/') {
        char drive = (char)std::tolower((unsigned char)win[0]);
        return "/mnt/" + std::string(1, drive) + win.substr(2);
    }

    // If already unix-y, just return
    return win;
}

// export CGAL coordinate x as an integer by multiplying by scale and rounding to long long 
// M2 wants integer coordinates 
// so we convert exact coordinates into large integers while approx. preserving geometry 
// ---------- scaled integer ----------
template <class FT>
static long long scaled_int(const FT& x, double scale) {
    const double xd = CGAL::to_double(x);
    return (long long) llround(xd * scale);
}


// run external command and capture its stdout so we can parse the true/ false result from M2
// ---------- capture stdout via popen/_popen ----------
static bool run_capture_stdout(const std::string& cmd, std::string& out) {
    #if defined(_WIN32)
        FILE* pipe = _popen(cmd.c_str(), "r");
    #else
        FILE* pipe = popen(cmd.c_str(), "r");
    #endif
        if (!pipe) return false;

        char buffer[4096];
        out.clear();
        while (fgets(buffer, sizeof(buffer), pipe)) out += buffer;

    #if defined(_WIN32)
        _pclose(pipe);
    #else
        pclose(pipe);
    #endif
        return true;
}

// parse M2 output: look for true/ false in output and set value accordingly
// M2 prints exactly one line: true or false, we do a substring search just in case 
static bool parse_m2_bool(const std::string& s, bool& value) {
    // M2 prints "true" or "false"
    if (s.find("true")  != std::string::npos) { value = true;  return true; }
    if (s.find("false") != std::string::npos) { value = false; return true; }
    return false;
}

// Data structure to export the CGAL triangulation into M2 script
// old_ids: list of vertex ids present in triangulation node
// old_to_to_new: we map from global ids to contiguous 0,...,k-1
// faces_new: triangles expressed using new contiguous indices
// ---------- export triangulation to (A, tri) with local indices ----------
struct ExportTri {
    std::vector<df::vertex_id> old_ids;                
    std::unordered_map<size_t, int> old_to_new;         
    std::vector<std::array<int,3>> faces_new;          
};

// interates over all finite vetices in a triangulation and builds 
// collects global vertex ids 
// sorts and makes them unique to get old_ids
// builds old_to_new map from global id to local index
// itertates over all finite faces 
// converts all vertex ids of faces into local indices and stores them 
static ExportTri extract_export_tri(const df::Tri2& T) {
    ExportTri out;

    std::vector<df::vertex_id> vids;
    for (auto v = T.finite_vertices_begin(); v != T.finite_vertices_end(); ++v) {
        vids.push_back(v->info());
    }
    std::sort(vids.begin(), vids.end());
    vids.erase(std::unique(vids.begin(), vids.end()), vids.end());

    out.old_ids = vids;
    out.old_to_new.reserve(vids.size() * 2);
    for (int i = 0; i < (int)vids.size(); ++i) {
        out.old_to_new[(size_t)vids[i]] = i;
    }

    for (auto f = T.finite_faces_begin(); f != T.finite_faces_end(); ++f) {
        const df::vertex_id a = f->vertex(0)->info();
        const df::vertex_id b = f->vertex(1)->info();
        const df::vertex_id c = f->vertex(2)->info();
        out.faces_new.push_back({
            out.old_to_new[(size_t)a],
            out.old_to_new[(size_t)b],
            out.old_to_new[(size_t)c]
        });
    }

    return out;
}

// writs a M2 file that defines point matrix A as integer matrix (columns are points) 
// definees triangultion tri as list of triangles {i,j,k} with local indices
// set A = triangulation (A, tri)
// sanity check if T is well defined in the sense of M2 triangulation 
// print result of regularity check 
// exit with 0
// ---------- write M2 script ----------
static std::string write_m2_script(
    const df::InputData& D,
    const ExportTri& X,
    double scale,
    const std::string& script_path
) {
    std::ofstream os(script_path, std::ios::binary);
    if (!os) return "failed to open script file for writing: " + script_path;

    os << "needsPackage \"Triangulations\";\n";

    // A = transpose matrix {{x0,y0},...}
    os << "A = transpose matrix {\n";
    for (int i = 0; i < (int)X.old_ids.size(); ++i) {
        const df::vertex_id old = X.old_ids[i];
        const auto& p = D.points2d[old]; 
        const long long xi = scaled_int(p.x(), scale);
        const long long yi = scaled_int(p.y(), scale);
        os << "  {" << xi << "," << yi << "}" << (i + 1 < (int)X.old_ids.size() ? "," : "") << "\n";
    }
    os << "};\n";

    os << "tri = {\n";
    for (int t = 0; t < (int)X.faces_new.size(); ++t) {
        const auto& f = X.faces_new[t];
        os << "  {" << f[0] << "," << f[1] << "," << f[2] << "}" << (t + 1 < (int)X.faces_new.size() ? "," : "") << "\n";
    }
    os << "};\n";

    os << "T = triangulation(A, tri);\n";
    os << "if (not isWellDefined T) then (print \"ERROR:notWellDefined\"; exit 2);\n";
    os << "print isRegularTriangulation T;\n";
    os << "exit 0;\n";

    return "";
}

// takes one node and checks regularity by writing M2 script and running it via WSL
static bool check_one_node_wsl_m2(
    const df::InputData& D,
    const df::Tri2& tri,
    double scale,
    const std::string& wsl_cmd,
    bool& is_regular,
    std::string& err) {
    // export triangulation into M2 script format
    ExportTri X = extract_export_tri(tri);

    // check if empty 
    if (X.faces_new.empty() || X.old_ids.size() < 3) {
        is_regular = true;
        err.clear();
        return true;
    }

    // write script into current working directory (Windows side path)
    const std::string script_win = tmp_filename("regcheck", ".m2");
    const std::string werr = write_m2_script(D, X, scale, script_win);
    if (!werr.empty()) { err = werr; return false; }

    // convert to WSL path
    // if script is in cwd, we can ask Windows for full path via relative, but simplest:
    // assume current working directory is shared (/mnt/c/...) and relative path works; still WSL runs in its own CWD
    // therefore we use absolute Windows path if possible
    // we use: /mnt/c + current dir is not known
    // write script into same folder where we launch .exe
    // if that folder is Windows path, WSL can access it by the translated path
    
   
    
    // resolve absolute path on Windows using _fullpath if available
    char absbuf[4096];
    #if defined(_WIN32)
        if (_fullpath(absbuf, script_win.c_str(), sizeof(absbuf)) == nullptr) {
            // fallback: relative
            std::snprintf(absbuf, sizeof(absbuf), "%s", script_win.c_str());
        }
    #else
        std::snprintf(absbuf, sizeof(absbuf), "%s", script_win.c_str());
    #endif
    const std::string script_wsl = win_to_wsl_path(std::string(absbuf));

    // run: wsl M2 -q <script_wsl>
    // quote the script path (spaces)
    std::ostringstream cmd;
    cmd << wsl_cmd << " M2 -q \"" << script_wsl << "\"";

    std::string out;
    if (!run_capture_stdout(cmd.str(), out)) {
        err = "failed to run command: " + cmd.str();
        return false;
    }

    // cleanup
    std::remove(script_win.c_str());

    if (out.find("ERROR:notWellDefined") != std::string::npos) {
        err = "M2: triangulation not well-defined (export issue / degeneracy / scale)";
        return false;
    }

    bool val = false;
    if (!parse_m2_bool(out, val)) {
        err = "M2 did not print true/false. Output:\n" + out;
        return false;
    }

    is_regular = val;
    err.clear();
    return true;
}

} // namespace


// final function to check regularity of all nodes in P1 
// prints only nodes which are not regular or where the check failed (with error message)
std::vector<NodeRegularity> check_poset_regularity_wsl(
    const df::InputData& D,
    const pst::Poset1& P1,
    double scale,
    const std::string& wsl_cmd,
    bool print_table) {

    std::vector<NodeRegularity> results;
    results.reserve(P1.nodes.size());

    int nonreg_count = 0;
    int fail_count = 0;

    for (int idx = 0; idx < (int)P1.nodes.size(); ++idx) {
        df::Tri2 tri = D.tri_upper;
        pst::replay_history_poset(tri, P1.nodes[idx].history, D);

        int numV = 0;
        for (auto v = tri.finite_vertices_begin(); v != tri.finite_vertices_end(); ++v) ++numV;

        int numF = 0;
        for (auto f = tri.finite_faces_begin(); f != tri.finite_faces_end(); ++f) ++numF;

        NodeRegularity r;
        r.node_idx = idx;
        r.num_vertices = numV;
        r.num_faces = numF;

        bool is_reg = false;
        std::string err;

        if (!check_one_node_wsl_m2(D, tri, scale, wsl_cmd, is_reg, err)) {
            r.error = err;
            ++fail_count;
        } else {
            r.is_regular = is_reg;
            if (!is_reg) ++nonreg_count;
        }

        results.push_back(std::move(r));
    }

    if (print_table) {
        std::cout << "\n=== Regularity check via WSL+Macaulay2 ===\n";
        std::cout << "nodes: " << results.size()
                  << " | non-regular: " << nonreg_count
                  << " | failed: " << fail_count << "\n";
        std::cout << "scale=" << scale << " | cmd=" << wsl_cmd << " M2\n\n";

        std::cout << std::left
                  << std::setw(8)  << "node"
                  << std::setw(10) << "|V|"
                  << std::setw(10) << "|F|"
                  << std::setw(12) << "regular?"
                  << "note\n";
        std::cout << "-------------------------------------------------------------\n";

        for (const auto& r : results) {
            if (r.error.empty() && r.is_regular) continue; 
            std::cout << std::left
                      << std::setw(8)  << r.node_idx
                      << std::setw(10) << r.num_vertices
                      << std::setw(10) << r.num_faces
                      << std::setw(12) << (r.error.empty() ? (r.is_regular ? "true" : "false") : "ERROR")
                      << (r.error.empty() ? (r.is_regular ? "" : "NON-REGULAR") : r.error)
                      << "\n";
        }
        std::cout << "-------------------------------------------------------------\n\n";
    }

    return results;
}


} // namespace regcheck
