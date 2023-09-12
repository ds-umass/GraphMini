//
// Created by ubuntu on 1/1/23.
//

#include "graph_converter.h"
#include "typedef.h"
#include "timer.h"
#include "meta.h"
#include "constant.h"
#include "logging.h"
#include <fstream>
#include <algorithm>
#include <oneapi/tbb/parallel_sort.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_reduce.h>
#include <oneapi/tbb/parallel_scan.h>

namespace minigraph {
    inline bool nextSNAPline(std::ifstream &infile, std::string &line, std::istringstream &iss,
                             uint64_t &src, uint64_t &dest) {
        do {
            if(!getline(infile, line)) return false;
        } while(line.length() == 0 || line[0] == '#');
        iss.clear();
        iss.str(line);
        return !!(iss >> src >> dest);
    }

    inline void getID(std::vector<uint64_t> &idMap, uint64_t &id, uint64_t &nextID) {
        if(idMap.size() <= id) {
            idMap.resize(id + Constant::kResizeIncrement, Constant::EmptyID<uint64_t>());
        }

        if(idMap.at(id) == Constant::EmptyID<uint64_t>()) {
            idMap.at(id) = nextID;
            nextID++;
        }
        id = idMap.at(id);
    }

    inline uint64_t parallel_max(const std::vector<uint64_t>& vec) {
        return tbb::parallel_reduce(
                tbb::blocked_range<uint64_t>(0, vec.size()),
                uint64_t(0),
                [&vec](const tbb::blocked_range<uint64_t>& r, uint64_t running_max){
                    for (uint64_t i = r.begin(); i < r.end(); ++i) {
                        running_max = std::max(vec[i], running_max);
                    }
                    return running_max;
                }, [](uint64_t a, uint64_t b){return std::max(a, b);} );
    }

    inline uint64_t parallel_sum(const std::vector<uint64_t>& vec) {
        return tbb::parallel_reduce(
                tbb::blocked_range<uint64_t>(0, vec.size()),
                uint64_t(0),
                [&vec](const tbb::blocked_range<uint64_t>& r, uint64_t running_total){
                    for (uint64_t i = r.begin(); i < r.end(); ++i) {
                        running_total += vec[i];
                    }
                    return running_total;
                }, std::plus<uint64_t>() );
    }

    inline std::vector<uint64_t> parallel_prefix_sum(const std::vector<uint64_t>& vec) {
        std::vector<uint64_t > res;
        res.resize(vec.size() + 1);
        tbb::parallel_scan(
                tbb::blocked_range<uint64_t>(0, vec.size()),
                        uint64_t (0),
                        [&vec, &res](const tbb::blocked_range<uint64_t>& r, uint64_t sum, bool is_final_scan) -> uint64_t
                {
                    uint64_t temp = sum;
                    for (uint64_t i = r.begin(); i < r.end(); i++){
                        temp += vec.at(i);
                        if (is_final_scan) res.at(i + 1) = temp;
                    }
                    return temp;
                },
                [] (uint64_t left, uint64_t right) {
                    return left + right;
                }
                );
        return res;
    }
    void GraphConverter::load_txt() {
        typedef std::pair<uint64_t , uint64_t> Edge;
        typedef std::vector<std::pair<uint64_t, uint64_t>> EdgeVector;
        LOG(INFO) << "Loading from file: " << data_file;

        Timer t;
        EdgeVector edge_vec;
        std::vector<uint64_t > idMap;
        uint64_t max_id = 0, src = 0, dst = 0, nextID = 0, line_num = 0;
        std::string line;
        std::istringstream iss;
        std::ifstream infile(data_file.c_str());
        while (nextSNAPline(infile, line, iss, src, dst)) {
            if (src == dst) continue;
            getID(idMap, src, nextID);
            getID(idMap, dst, nextID);
            max_id = std::max(max_id, src);
            max_id = std::max(max_id, dst);
            if (max_id >= degrees.size()) degrees.resize(max_id + Constant::kResizeIncrement, 0);
            edge_vec.push_back(std::make_pair(src, dst)); // make it undirected
            edge_vec.push_back(std::make_pair(dst, src));
            line_num++;
            if (line_num % int(1e6) == 0){
                LOG(INFO) << "Read: " << line_num / int(1e6) << "M edges";
            }
        }

        LOG(INFO) << "Finished Reading in: " << t.Passed() << " seconds";
        // make sort edge pair for building indices
        t.Reset();
        tbb::parallel_sort(edge_vec.begin(), edge_vec.end(),
                  [](const Edge& a, const Edge& b)->bool {
                      if (a.first == b.first) {
                          return a.second < b.second;
                      } else {
                          return a.first < b.first;
                      }
                  });
        edge_vec.erase(std::unique(edge_vec.begin(), edge_vec.end()), edge_vec.end());
        edge_vec.shrink_to_fit();

        LOG(INFO) << "Finished sorting edges: " << t.Passed() << " seconds";
        t.Reset();

        v_num = nextID;
        // build indices
        degrees.resize(v_num, 0);
        offsets.resize(v_num, 0);
        triangles.resize(v_num, 0);

        indices.resize(edge_vec.size(), Constant::EmptyID<uint64_t>());
        for (uint64_t i = 0; i < edge_vec.size(); ++i) {
            const auto& pair = edge_vec.at(i);
            indices.at(i) = pair.second;
            degrees.at(pair.first)++;
        }
        indptr = parallel_prefix_sum(degrees);
        edge_vec.clear();
        edge_vec.shrink_to_fit();
        e_num = indices.size();


        LOG(INFO) << "Finished building indices in: " << t.Passed() << " seconds";
        tbb::parallel_for(tbb::blocked_range<uint64_t >(0, v_num), [this](tbb::blocked_range<uint64_t> r){
            for (uint64_t i = r.begin(); i != r.end(); i++){
                Counter counter{};
                auto v1_start = indices.cbegin() + indptr.at(i);
                auto v1_end = indices.cbegin() + indptr.at(i+1);
                offsets.at(i) = std::distance(v1_start, std::lower_bound(v1_start, v1_end, i));

                const int v1_degree = std::distance(v1_start, v1_end);
                for (int j = 0; j < v1_degree; j++){
                    uint64_t v2_id = v1_start[j];
                    auto v2_start = indices.cbegin() + indptr.at(v2_id);
                    auto v2_end = indices.cbegin() + indptr.at(v2_id+1);
                    std::set_intersection(v1_start, v1_end, v2_start, v2_end, std::back_inserter(counter));
                }
                triangles.at(i) = counter.count;
            }
        });
        tri_num = parallel_sum(triangles);
        max_deg = parallel_max(degrees);
        max_offset = parallel_max(offsets);
        max_tri = parallel_max(triangles);
        LOG(INFO) << "Finished counting triangles in: " << t.Passed() << " seconds";
    }


    void GraphConverter::save_meta() {
        MetaData meta(v_num, e_num, tri_num, max_deg, max_offset, max_tri);
        meta.save(in_dir);
    }

    std::vector<uint32_t> to_32(const std::vector<uint64_t>& data) {
        std::vector<uint32_t > out;
        out.resize(data.size());
        for (uint64_t i = 0; i < data.size(); ++i) {
            out.at(i) = static_cast<uint32_t>(data.at(i));
        }
        return out;
    }

    void save_u32(std::filesystem::path path, const std::vector<uint32_t>& data) {
        std::ofstream out;
        out.open(path, std::ios::binary | std::ios::out);
        out.write(reinterpret_cast<const char *>(&data[0]), sizeof(uint32_t) * data.size() );
        out.close();
    }

    void save_u64(std::filesystem::path path, const std::vector<uint64_t>& data) {
        std::ofstream out;
        out.open(path, std::ios::binary | std::ios::out);
        out.write(reinterpret_cast<const char *>(&data[0]), sizeof(uint64_t) * data.size() );
        out.close();
    }

    void GraphConverter::save_bin_u32() {
        std::filesystem::path indicesPath = in_dir / Constant::kIndicesU32File;
        std::filesystem::path indptrPath = in_dir / Constant::kIndptrU32File;
        std::filesystem::path trianglePath = in_dir / Constant::kTriangleU32File;
        std::filesystem::path degreePath = in_dir / Constant::kDegreeU32File;
        std::filesystem::path offsetPath = in_dir / Constant::kOffsetU32File;
        if (v_num < std::numeric_limits<uint32_t>::max()) save_u32(indicesPath, to_32(indices));
        if (e_num < std::numeric_limits<uint32_t>::max()) save_u32(indptrPath, to_32(indptr));
        if (max_deg < std::numeric_limits<uint32_t>::max()) save_u32(degreePath, to_32(degrees));
        if (max_offset < std::numeric_limits<uint32_t>::max()) save_u32(offsetPath, to_32(offsets));
        if (max_tri < std::numeric_limits<uint32_t>::max()) save_u32(trianglePath, to_32(triangles));
    }

    void GraphConverter::save_bin() {
        Timer t;
        std::filesystem::path indicesPath = in_dir / Constant::kIndicesU64File;
        std::filesystem::path indptrPath = in_dir / Constant::kIndptrU64File;
        std::filesystem::path trianglePath = in_dir / Constant::kTriangleU64File;
        std::filesystem::path degreePath = in_dir / Constant::kDegreeU64File;
        std::filesystem::path offsetPath = in_dir / Constant::kOffsetU64File;
        save_u64(indicesPath, indices);
        save_u64(indptrPath, indptr);
        save_u64(trianglePath, triangles);
        save_u64(degreePath, degrees);
        save_u64(offsetPath, offsets);

//        std::ofstream outfile;
//        outfile.open(indicesPath, std::ios::binary | std::ios::out);
//        outfile.write( reinterpret_cast<const char *>(&indices[0]), sizeof(uint64_t) * indices.size());
//        outfile.close();
//
//        outfile.open(indptrPath, std::ios::binary | std::ios::out);
//        outfile.write( reinterpret_cast<const char *>(&indptr[0]), sizeof(uint64_t) * indptr.size());
//        outfile.close();
//
//        outfile.open(trianglePath, std::ios::binary | std::ios::out);
//        outfile.write( reinterpret_cast<const char *>(&triangles[0]), sizeof(uint64_t) * triangles.size());
//        outfile.close();
//
//        outfile.open(degreePath, std::ios::binary | std::ios::out);
//        outfile.write( reinterpret_cast<const char *>(&degrees[0]), sizeof(uint64_t) * degrees.size());
//        outfile.close();
//
//        outfile.open(offsetPath, std::ios::binary | std::ios::out);
//        outfile.write( reinterpret_cast<const char *>(&offsets[0]), sizeof(uint64_t) * offsets.size());
//        outfile.close();

        save_bin_u32();
        LOG(INFO) << "Finished writing binary files in: " << t.Passed() << " seconds";
    }

    void GraphConverter::convert(std::filesystem::path input_dir) {
        in_dir = input_dir;
        data_file = in_dir / "snap.txt";
        CHECK(std::filesystem::is_regular_file(data_file)) << "Cannot find file: " << data_file;

        indices.clear();
        indptr.clear();
        triangles.clear();
        offsets.clear();
        degrees.clear();

        load_txt();
        save_meta();
        save_bin();
    }
}
