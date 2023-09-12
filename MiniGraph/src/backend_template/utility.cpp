//
// Created by ubuntu on 1/19/23.
//

#include "utility.h"
#include "meta.h"
#include <filesystem>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>
#include <iostream>
namespace minigraph
{
    template<typename T>
    void read_file(std::filesystem::path path, T *& pointer, uint64_t num_elements){
        CHECK(std::filesystem::is_regular_file(path)) << "File does not exists: " << path;
        if (pointer != nullptr) free(pointer);
        const size_t num_bytes = sizeof(T) * num_elements;
        pointer = new T[num_elements];
        std::ifstream file;
        file.open(path, std::ios::binary | std::ios::in);
        file.read(reinterpret_cast<char *>(pointer), num_bytes);
        CHECK(file.gcount() == (long int) num_bytes) << "Only read " << ToReadableSize(file.gcount()) << " out of "<< ToReadableSize(num_bytes) <<" from " << path;
        file.close();
        LOG(INFO) << "Read " << ToReadableSize(num_bytes) << " from " << path;
    };

    template<typename T>
    void mmap_file(std::filesystem::path path, T *& pointer, uint64_t num_elements){
        CHECK(std::filesystem::is_regular_file(path)) << "File does not exists: " << path;
        if (pointer != nullptr) free(pointer);
        int fd = open(path.c_str(), O_RDONLY, 0);
        CHECK(fd != -1) << "Failed to open: " << path;
        pointer = static_cast<T*>(
                mmap(nullptr, sizeof(T) * num_elements, PROT_READ, MAP_SHARED, fd, 0));
        CHECK(pointer != MAP_FAILED) << "Failed to map file: " << path;
        CHECK(close(fd) == 0) << "Failed to close file: " << path;
    };

    template<class IdType>
    Graph<IdType>* load_bin(std::string _in_dir, bool _mmap) {
        Graph<IdType>* out = new Graph<IdType>;
        MetaData m_meta;
        m_meta.read(_in_dir);

        out->m_mmap = _mmap;
        out->num_vertex = m_meta.num_vertex;
        out->num_edge = m_meta.num_edge;
        out->num_triangle = m_meta.num_triangle;
        out->max_offset = m_meta.max_offset;
        out->max_degree = m_meta.max_degree;
        out->max_triangle = m_meta.max_triangle;

        std::filesystem::path indicesFile = _in_dir;
        if (sizeof(IdType) == sizeof(uint64_t)) indicesFile /= Constant::kIndicesU64File;
        else if(sizeof(IdType) == sizeof(uint32_t)) indicesFile /= Constant::kIndicesU32File;
        else exit(-1 && "unsupported IdType");

        read_file<uint64_t>(std::filesystem::path{_in_dir} / Constant::kIndptrU64File, out->m_indptr, m_meta.num_vertex + 1);
        read_file<uint64_t>(std::filesystem::path{_in_dir} / Constant::kOffsetU64File, out->m_offset, m_meta.num_vertex);
        read_file<uint64_t>(std::filesystem::path{_in_dir} / Constant::kTriangleU64File, out->m_triangles, m_meta.num_vertex);
        if (_mmap) {
            mmap_file<IdType>(indicesFile, out->m_indices, m_meta.num_edge);
        } else {
            read_file<IdType>(indicesFile, out->m_indices, m_meta.num_edge);
        }
        return out;
    }
    template Graph<uint32_t>* load_bin<uint32_t>(std::string _in_dir, bool _mmap);
    template Graph<uint64_t>* load_bin<uint64_t>(std::string _in_dir, bool _mmap);

    template<typename IdType>
    std::ostream &operator<<(std::ostream &os, const VertexSet<IdType> &dt) {
        if (dt.vid() == Constant::EmptyID<IdType>()) {
            os << "VertexSet(-1)\t=\t[";
        } else {
            os << "VertexSet(" << dt.vid() << ")\t=\t[";
        }
        constexpr uint64_t expected_to_show = 10;
        uint64_t actual_show = std::min(expected_to_show, dt.size());
        for(size_t i = 0; i < actual_show; i++){
            os << dt[i] << ",";
        }
        if (expected_to_show < dt.size()) os << "...";
        os << "]\t" ;
        os << "(size=" << dt.size() << " pooled="<< dt.pooled() << ")\n";
        return os;
    }

    template std::ostream &operator<< <uint32_t> (std::ostream &os, const VertexSet<uint32_t> &dt);
    template std::ostream &operator<< <uint64_t> (std::ostream &os, const VertexSet<uint64_t> &dt);
}