//
// Created by ubuntu on 2/2/23.
//

#ifndef MINIGRAPH_GRAPH_H
#define MINIGRAPH_GRAPH_H
#include "vertex_set.h"
#include <sys/mman.h>
namespace minigraph {
    struct Graph {
        IdType *m_indices{nullptr};
        uint64_t *m_indptr{nullptr};
        uint64_t *m_offset{nullptr};
        uint64_t *m_triangles{nullptr};
        uint64_t num_vertex{0}, num_edge{0}, num_triangle{0};
        uint64_t max_degree{0}, max_offset{0}, max_triangle{0};
        bool m_mmap{false};

        Graph() = default;

        ~Graph() {
            if (!m_mmap) {
                if (m_indices != nullptr) delete[] m_indices;
            } else {
                munmap(m_indices, sizeof(IdType) * num_edge);
            }
            if (m_indptr != nullptr) delete[] m_indptr;
            if (m_offset != nullptr) delete[] m_offset;
            if (m_triangles != nullptr) delete[] m_triangles;
        };

        uint64_t get_vnum() const { return num_vertex; }

        uint64_t get_enum() const { return num_edge; };

        uint64_t get_tnum() const { return num_triangle; };

        uint64_t get_maxdeg() const { return max_triangle; }

        uint64_t get_maxoffset() const { return max_offset; }

        uint64_t get_maxtri() const { return max_triangle; }

        uint64_t Degree(IdType v_id) const {
            assert(v_id < num_vertex);
            return m_indptr[v_id + 1] - m_indptr[v_id]; };

        uint64_t Offset(IdType v_id) const { assert(v_id < num_vertex); return m_offset[v_id]; };

        // return adj of v
        VertexSet N(IdType v_id) const {
            auto start = m_indices + m_indptr[v_id];
            auto degree = Degree(v_id);
            return VertexSet(v_id, start, degree);
        };

        // return adj of v bounded by v_id
        VertexSet NBound(IdType v_id) const {
            auto start = m_indices + m_indptr[v_id];
            auto degree = Offset(v_id);
            return VertexSet(v_id, start, degree);
        };
    };


}
#endif //MINIGRAPH_GRAPH_H
