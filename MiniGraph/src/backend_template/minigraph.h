//
// Created by ubuntu on 1/1/23.
//
#ifndef MINIGRAPH_MINIGRAPH_H
#define MINIGRAPH_MINIGRAPH_H
#include "common.h"
#include "graph.h"
#include "vertex_set.h"
#include <vector>
#include <oneapi/tbb/parallel_for.h>
#include <cstring>
#include <cstdlib>
namespace minigraph {
    using IdType = uint32_t;
    constexpr IdType INVALID_ID = static_cast<IdType>(-1);
    struct Container {
        IdType* data{nullptr};
        size_t capacity{0}; // number of elements can be stored into the buffer
        Container() = default;
        Container(IdType *_data, size_t _capacity):
                data{_data}, capacity{_capacity}{};

        IdType& operator[](size_t i) {
            assert(i < capacity);
            return data[i];
        };

        const IdType& operator[](size_t i) const {
            assert(i < capacity);
            return data[i];
        }
    };

    class MiniGraphPool
    {
    private:
        std::vector<Container> buffer_exist; // large buffer (>4kb)
        std::vector<Container> buffer_avail; // large buffer (>4kb)

        std::vector<Container> m_4k_buffer_exist; // 4kb pages buffer
        std::vector<Container> m_4k_buffer_avail; // 4kb pages buffer
        constexpr static int m_increment = 4096 / sizeof(IdType);
        constexpr static int m_4k_cap = 4096 / sizeof(IdType);

        Container AllocateSmallWorkSpace(){
            if (m_4k_buffer_avail.empty()) {
                IdType * data = new IdType [m_4k_cap];
                Container out{data, m_4k_cap};
                buffer_exist.push_back(out);
                TOTAL_ALLOCATED += m_4k_cap * sizeof(IdType);
                return out;
            } else {
                Container out = m_4k_buffer_avail.back();
                m_4k_buffer_avail.pop_back();
                return out;
            }
        }

        void FreeSmallWorkSpace(Container ctn) {
            m_4k_buffer_avail.push_back(ctn);
        }
    public:
        inline static std::atomic_uint64_t TOTAL_ALLOCATED{0};
        MiniGraphPool() = default;
        ~MiniGraphPool() {
            for (auto x: buffer_exist) {
                delete[] x.data;
            }

            for (auto x: m_4k_buffer_exist) {
                delete[] x.data;
            }
        }

        static MiniGraphPool& Get() {
            static thread_local MiniGraphPool pool;
            return pool;
        }

        Container AllocateWorkSpace(size_t num_elements){
            if (num_elements <= m_4k_cap) return AllocateSmallWorkSpace();
            if (buffer_avail.empty() || buffer_avail.back().capacity < num_elements)
            {
                size_t capacity = ((num_elements / m_increment) + 1) * m_increment;
                IdType * data = new IdType[capacity];
                Container ctn(data, capacity);
                buffer_exist.push_back(ctn);
                TOTAL_ALLOCATED += capacity * sizeof(IdType);
                return ctn;
            }
            for (auto itr = buffer_avail.begin(); itr!= buffer_avail.end(); itr++){
                if (itr->capacity >= num_elements){
                    Container out = *itr;
                    buffer_avail.erase(itr);
                    return out;
                }
            }
            LOG(FATAL) << "Logic error (should not reach here)";
            exit(-1);
        }

        void FreeWorkSpace(Container ctn){
            if (ctn.data == nullptr) return;
            if (ctn.capacity <= m_4k_cap) return FreeSmallWorkSpace(ctn);
            for (auto itr = buffer_avail.begin(); itr != buffer_avail.end(); itr++) {
                if(itr->capacity >= ctn.capacity){
                    buffer_avail.insert(itr, ctn);
                    return;
                }
            }
            buffer_avail.push_back(ctn);
        }

        void Resize(Container& ctn, size_t new_capacity){
            if (ctn.data == nullptr) ctn = AllocateWorkSpace(new_capacity);
            Container out = AllocateWorkSpace(new_capacity);
            std::copy(ctn.data, ctn.data + ctn.capacity, out.data);
            FreeWorkSpace(ctn);
            ctn = out;
        }
    };

    class ManagedContainer {
    private:
        Container m_ctn{nullptr, 0};
        size_t m_size{0};
    public:
        ManagedContainer() = default;
        explicit ManagedContainer(size_t _size) {
            m_ctn = MiniGraphPool::Get().AllocateWorkSpace(_size);
            m_size = _size;
        }

        ~ManagedContainer() {
            MiniGraphPool::Get().FreeWorkSpace(m_ctn);
        }

        ManagedContainer(const ManagedContainer& src) = delete;
        ManagedContainer &operator=(const ManagedContainer& src) = delete;
        void swap(ManagedContainer& src) noexcept {
            std::swap(m_ctn, src.m_ctn);
            std::swap(m_size, src.m_size);
        }
        ManagedContainer &operator=(ManagedContainer &&src) noexcept {
            swap(src); return *this;
        }
        ManagedContainer (ManagedContainer &&src) noexcept { swap(src);};
        void set_size(size_t _size) {m_size = _size;};
        size_t size() const {return m_size;};
        size_t capacity() const {return m_ctn.capacity;};
        IdType *begin(){ return m_ctn.data; }
        const IdType *begin() const { return m_ctn.data; }
        IdType *end() { return m_ctn.data + m_size; }
        const IdType *end() const { return m_ctn.data + m_size; }
        IdType & operator[] (size_t i) {return m_ctn[i];}
        const IdType & operator[] (size_t i) const {return m_ctn[i];}
        void Resize(size_t _capacity){
            MiniGraphPool::Get().Resize(m_ctn, _capacity);
        };
        void Reserve(size_t _capacity){
            MiniGraphPool::Get().FreeWorkSpace(m_ctn);
            m_ctn = MiniGraphPool::Get().AllocateWorkSpace(_capacity);
        };
    };

    inline ManagedContainer get_indices(const VertexSet<IdType>& _vertices, const VertexSet<IdType>& _to_iter){
        ManagedContainer out(_to_iter.size());
        size_t idx_l = 0, idx_r = 0, idx_out = 0;
        size_t _v_size = _vertices.size();
        size_t _iter_size = _to_iter.size();
        while(idx_l < _v_size && idx_r < _iter_size) {
            const IdType left = _vertices[idx_l];
            const IdType right = _to_iter[idx_r];
            if(left == right) out[idx_out++] = idx_l;
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
        }
        assert(out.size() == _to_iter.size() && "Logic error");
        out.set_size(idx_out);
        return out;
    };

    struct MiniGraphIF
    {
        MiniGraphIF() = default;
        virtual ~MiniGraphIF() = default;
        inline static const Graph<IdType> * DATA_GRAPH{nullptr};

        virtual void build(const VertexSet<IdType>& _vertex,
                           const VertexSet<IdType>& _intersect,
                           const VertexSet<IdType>& _iter) = 0;

        virtual void build(MiniGraphIF& mg,
                           const VertexSet<IdType>& _vertex,
                           const VertexSet<IdType>& _intersect,
                           const VertexSet<IdType>& _iter) = 0;

        virtual VertexSet<IdType> N(IdType i) = 0;
        virtual ManagedContainer indices(const VertexSet<IdType>& _to_iter) const = 0;
        virtual IdType Degree(IdType i) const = 0;
    };

    class MiniGraphEager : public MiniGraphIF
    {
    private:
        VertexSet<IdType> m_vertex;
        VertexSet<IdType> m_intersect;
        ManagedContainer m_ctn;
        ManagedContainer m_pos;
        ManagedContainer m_degree;
        const bool m_bounded{false};
        const bool m_par{false};
        size_t num_edges{0};
    public:
        MiniGraphEager() = default;
        ~MiniGraphEager() = default;
        MiniGraphEager(bool _bounded, bool _par): m_bounded{_bounded}, m_par{_par}{};
        void build(const VertexSet<IdType>& _vertex,
                   const VertexSet<IdType>& _intersect,
                   const VertexSet<IdType>& _iter) override {
            m_vertex = _vertex;
            m_intersect = _intersect;
            if (m_pos.capacity() <= m_vertex.size()){
                m_pos.Reserve(m_vertex.size() + 1024);
                m_degree.Reserve(m_vertex.size() + 1024);
            }
            m_pos.set_size(m_vertex.size()+1);
            m_degree.set_size(m_vertex.size());
            num_edges = 0;
            m_pos[0] = num_edges;

            if(m_par){
                for (size_t i = 0; i < m_vertex.size(); i++){
                    IdType v_id = m_vertex[i];
                    m_pos[i + 1] = m_bounded ? m_pos[i] + std::min(DATA_GRAPH->Offset(v_id), m_intersect.size())
                                             : m_pos[i] + std::min(DATA_GRAPH->Degree(v_id), m_intersect.size());
                }
                IdType est_edges = m_pos[m_vertex.size()];
                if (m_ctn.capacity() < est_edges) m_ctn.Reserve(est_edges);

                tbb::parallel_for(tbb::blocked_range<uint64_t >(0, m_vertex.size()), [this](tbb::blocked_range<uint64_t> r){
                    for (uint64_t i = r.begin(); i < r.end(); ++i) {
                        IdType v_id = m_vertex[i];
                        IdType* start = m_ctn.begin() + m_pos[i];
                        size_t degree = m_bounded ? m_intersect.intersect(DATA_GRAPH->NBound(v_id), start)
                                                : m_intersect.intersect(DATA_GRAPH->N(v_id), start);
                        m_degree[i] = degree;
                    }
                });
            } else {
                for (uint64_t i = 0; i < m_vertex.size(); ++i) {
                    size_t buffer_required = m_intersect.size() + m_pos[i];
                    if (m_ctn.capacity() < buffer_required) m_ctn.Resize(buffer_required);
                    IdType v_id = m_vertex[i];
                    IdType* start = m_ctn.begin() + m_pos[i];
                    size_t degree = m_bounded ? m_intersect.intersect(DATA_GRAPH->NBound(v_id), start)
                                            : m_intersect.intersect(DATA_GRAPH->N(v_id), start);
                    m_degree[i] = degree;
                    m_pos[i + 1] = m_pos[i] + degree;
                    num_edges+=degree;
                }
            }
        }

        void build(MiniGraphIF& mg,
                   const VertexSet<IdType>& _vertex,
                   const VertexSet<IdType>& _intersect,
                   const VertexSet<IdType>& _iter) override {
            m_vertex = _vertex;
            m_intersect = _intersect;
            if (m_pos.capacity() <= m_vertex.size()){
                m_pos.Reserve(m_vertex.size() + 1024);
                m_degree.Reserve(m_vertex.size() + 1024);
            }
            m_pos.set_size(m_vertex.size()+1);
            m_degree.set_size(m_vertex.size());
            m_pos[0] = 0;
            ManagedContainer m_indices = mg.indices(m_vertex);
            if(m_par){
                for (size_t i = 0; i < m_vertex.size(); i++){
                    IdType adj_idx = m_indices[i];
                    m_pos[i + 1] = m_pos[i] + std::min(mg.Degree(adj_idx), (IdType) m_intersect.size());
                }
                IdType est_edges = m_pos[m_vertex.size()];
                if (m_ctn.capacity() < est_edges) m_ctn.Reserve(est_edges);

                tbb::parallel_for(tbb::blocked_range<uint64_t >(0, m_vertex.size()), [this, &m_indices, &mg](tbb::blocked_range<uint64_t> r){
                    for (uint64_t i = r.begin(); i < r.end(); ++i) {
                        IdType v_id = m_vertex[i];
                        IdType adj_idx = m_indices[i];
                        IdType* start = m_ctn.begin() + m_pos[i];
                        size_t degree = m_bounded ? m_intersect.intersect(mg.N(adj_idx), v_id, start)
                                                  : m_intersect.intersect(mg.N(adj_idx), start);
                        m_degree[i] = degree;
                    }
                });
            } else {
                for (uint64_t i = 0; i < m_vertex.size(); ++i) {
                    size_t buffer_required = m_intersect.size() + m_pos[i];
                    if (m_ctn.capacity() < buffer_required) m_ctn.Resize(buffer_required);
                    IdType v_id = m_vertex[i];
                    IdType adj_idx = m_indices[i];
                    IdType* start = m_ctn.begin() + m_pos[i];
                    size_t degree = m_bounded ? m_intersect.intersect(mg.N(adj_idx), v_id, start)
                                              : m_intersect.intersect(mg.N(adj_idx), start);
                    m_degree[i] = degree;
                    m_pos[i + 1] = m_pos[i] + degree;
                }
            }
        }

        VertexSet<IdType> N(IdType i) override {
            return VertexSet<IdType>(Constant::EmptyID<IdType>(), m_ctn.begin() + m_pos[i], m_degree[i]);
        }

        ManagedContainer indices(const VertexSet<IdType>& _to_iter) const override {
            return get_indices(m_vertex, _to_iter);
        }
        IdType Degree(IdType i) const override {
            return m_degree[i];
        }
    };

    class MiniGraphLazy: public MiniGraphIF {
    private:
        VertexSet<IdType> m_vertex;
        VertexSet<IdType> m_intersect;
        ManagedContainer m_ctn;
        ManagedContainer m_pos;
        ManagedContainer m_degree;
        ManagedContainer m_indices;
        ManagedContainer m_mg_indices;
        MiniGraphIF* m_mg{nullptr};
        const bool m_bounded{false};
        const bool m_par{false};
        size_t num_edges{0};

    public:
        MiniGraphLazy() = default;
        MiniGraphLazy(bool _bounded, bool _par) : m_bounded{_bounded}, m_par{_par}{};
        ~MiniGraphLazy() = default;
        void build(const VertexSet<IdType> &_vertex,
                   const VertexSet<IdType> &_intersect,
                   const VertexSet<IdType> &_iter) override {
            m_vertex = _vertex;
            m_intersect = _intersect;

            if (m_pos.capacity() <= m_vertex.size()){
                m_pos.Reserve(m_vertex.size() + 1024);
                m_degree.Reserve(m_vertex.size() + 1024);
            }
            m_pos.set_size(m_vertex.size()+1);
            m_degree.set_size(m_vertex.size());
            num_edges = 0;
            m_pos[0] = 0;

            for (auto &x: m_degree) {
                x = INVALID_ID;
            }
            m_indices = get_indices(m_vertex, _iter);

            if(m_par){
                for (auto v_idx : m_indices){
                    IdType v_id = m_vertex[v_idx];
                    m_pos[v_idx] = num_edges;
                    num_edges += m_bounded ? std::min(DATA_GRAPH->Offset(v_id), m_intersect.size())
                                           : std::min(DATA_GRAPH->Degree(v_id), m_intersect.size());
                }
                if (m_ctn.capacity() < num_edges) m_ctn.Reserve(num_edges);

                tbb::parallel_for(tbb::blocked_range<uint64_t >(0, m_indices.size()), [this](tbb::blocked_range<uint64_t> r){
                    for (uint64_t i = r.begin(); i < r.end(); ++i) {
                        IdType v_idx = m_indices[i];
                        IdType v_id = m_vertex[v_idx];
                        IdType* start = m_ctn.begin() + m_pos[i];
                        size_t degree = m_bounded ? m_intersect.intersect(DATA_GRAPH->NBound(v_id), start)
                                                  : m_intersect.intersect(DATA_GRAPH->N(v_id), start);
                        m_degree[v_idx] = degree;
                    }
                });
            } else {
                for (auto v_idx : m_indices) {
                    m_pos[v_idx] = num_edges;
                    size_t buffer_required = m_intersect.size() + num_edges;
                    if (m_ctn.capacity() < buffer_required) m_ctn.Resize(buffer_required);
                    IdType v_id = m_vertex[v_idx];
                    IdType* start = m_ctn.begin() + m_pos[v_idx];
                    size_t degree = m_bounded ? m_intersect.intersect(DATA_GRAPH->NBound(v_id), start)
                                              : m_intersect.intersect(DATA_GRAPH->N(v_id), start);
                    m_degree[v_idx] = degree;
                    num_edges += degree;
                }
            }
        }

        void build(MiniGraphIF& mg,
                   const VertexSet<IdType> &_vertex,
                   const VertexSet<IdType> &_intersect,
                   const VertexSet<IdType> &_iter) override {
            m_vertex = _vertex;
            m_intersect = _intersect;
            m_mg = &mg;
            if (m_pos.capacity() <= m_vertex.size()){
                m_pos.Reserve(m_vertex.size() + 1024);
                m_degree.Reserve(m_vertex.size() + 1024);
            }
            m_pos.set_size(m_vertex.size()+1);
            m_degree.set_size(m_vertex.size());
            num_edges = 0;
            m_pos[0] = 0;
            for (auto &x: m_degree) {
                x = INVALID_ID;
            }
            m_indices = get_indices(m_vertex, _iter);
            m_mg_indices = mg.indices(m_vertex);
            assert(m_indices.size() == _iter.size());
            assert(m_mg_indices.size() == m_vertex.size());

            if(m_par){
                for (IdType v_idx : m_indices){
                    IdType adj_idx = m_mg_indices[v_idx];
                    m_pos[v_idx] = num_edges;
                    num_edges += std::min(mg.Degree(adj_idx), (IdType) m_intersect.size());
                }
                if (m_ctn.capacity() < num_edges) m_ctn.Reserve(num_edges);

                tbb::parallel_for(tbb::blocked_range<uint64_t >(0, m_indices.size()), [this](tbb::blocked_range<uint64_t> r){
                    for (uint64_t i = r.begin(); i < r.end(); ++i) {
                        IdType v_idx = m_indices[i];
                        IdType v_id = m_vertex[v_idx];
                        IdType adj_idx = m_mg_indices[v_idx];
                        IdType* start = m_ctn.begin() + m_pos[i];
                        size_t degree = m_bounded ? m_intersect.intersect(m_mg->N(adj_idx), v_id, start)
                                                  : m_intersect.intersect(m_mg->N(adj_idx), start);
                        m_degree[v_idx] = degree;
                    }
                });
            } else {
                for (IdType v_idx: m_indices) {
                    IdType v_id = m_vertex[v_idx];
                    IdType adj_idx = m_mg_indices[v_idx];
                    size_t buffer_required = m_intersect.size() + num_edges;
                    if (m_ctn.capacity() < buffer_required) m_ctn.Resize(buffer_required);
                    m_pos[v_idx] = num_edges;
                    IdType* start = m_ctn.begin() + m_pos[v_idx];
                    size_t degree = m_bounded ? m_intersect.intersect(m_mg->N(adj_idx), v_id, start)
                                              : m_intersect.intersect(m_mg->N(adj_idx), start);
                    m_degree[v_idx] = degree;
                    num_edges += degree;
                }
            }
        }

        VertexSet<IdType> N(IdType i) override {
            if (m_degree[i] == INVALID_ID) {
                if (m_mg) return m_mg->N(m_mg_indices[i]);
                else return DATA_GRAPH->N(m_vertex[i]);
            } else {
                return VertexSet<IdType>(Constant::EmptyID<IdType>(), m_ctn.begin() + m_pos[i], m_degree[i]);
            }
        }

        IdType Degree(IdType i) const override {
            return m_degree[i];
        }

        ManagedContainer indices(const VertexSet<IdType>& _to_iter) const override {
            return get_indices(m_vertex, _to_iter);
        }
    };



    class MiniGraphOnline: public MiniGraphIF {
    private:
        VertexSet<IdType> m_vertex;
        VertexSet<IdType> m_intersect;
        ManagedContainer m_ctn;
        ManagedContainer m_pos;
        ManagedContainer m_degree;
        ManagedContainer m_indices;
        ManagedContainer m_mg_indices;
        MiniGraphIF* m_mg{nullptr};
        const bool m_bounded{false};
        const bool m_par{false};
        size_t num_edges{0};
        size_t est_edges{0};
    public:
        MiniGraphOnline() = default;
        MiniGraphOnline(bool _bounded, bool _par) : m_bounded{_bounded}, m_par{_par}{};
        ~MiniGraphOnline() = default;
        void build(const VertexSet<IdType> &_vertex,
                   const VertexSet<IdType> &_intersect,
                   const VertexSet<IdType> &_iter) override {
            m_vertex = _vertex;
            m_intersect = _intersect;

            if (m_pos.capacity() <= m_vertex.size()){
                m_pos.Reserve(m_vertex.size() + 1024);
                m_degree.Reserve(m_vertex.size() + 1024);
            }
            m_pos.set_size(m_vertex.size()+1);
            m_degree.set_size(m_vertex.size());
            num_edges = 0;
            m_pos[0] = 0;
            for (auto &x: m_degree) {
                x = INVALID_ID;
            }
            m_indices = get_indices(m_vertex, _iter);
            for (size_t i = 0; i < m_vertex.size(); i++){
                IdType v_id = m_vertex[i];
                m_pos[i + 1] = m_bounded ? m_pos[i] + std::min(DATA_GRAPH->Offset(v_id), m_intersect.size())
                                         : m_pos[i] + std::min(DATA_GRAPH->Degree(v_id), m_intersect.size());
            }
            est_edges = m_pos[m_vertex.size()];
            if (m_ctn.capacity() < est_edges) m_ctn.Reserve(est_edges);
            if (m_par){
                tbb::parallel_for(tbb::blocked_range<uint64_t >(0, m_indices.size()), [this](tbb::blocked_range<uint64_t> r){
                    for (uint64_t i = r.begin(); i < r.end(); ++i) {
                        IdType v_idx = m_indices[i];
                        IdType v_id = m_vertex[v_idx];
                        IdType* start = m_ctn.begin() + m_pos[v_idx];
                        size_t degree = m_bounded ? m_intersect.intersect(DATA_GRAPH->NBound(v_id), start)
                                                  : m_intersect.intersect(DATA_GRAPH->N(v_id), start);
                        m_degree[v_idx] = degree;
                    }
                });
            } else {
                for (auto v_idx : m_indices) {
                    IdType v_id = m_vertex[v_idx];
                    IdType* start = m_ctn.begin() + m_pos[v_idx];
                    size_t degree = m_bounded ? m_intersect.intersect(DATA_GRAPH->NBound(v_id), start)
                                              : m_intersect.intersect(DATA_GRAPH->N(v_id), start);
                    m_degree[v_idx] = degree;
                    num_edges += degree;
                }
            }
        }

        void build(MiniGraphIF& mg,
                   const VertexSet<IdType> &_vertex,
                   const VertexSet<IdType> &_intersect,
                   const VertexSet<IdType> &_iter) override {
            m_mg = &mg;
            m_vertex = _vertex;
            m_intersect = _intersect;

            if (m_pos.capacity() <= m_vertex.size()){
                m_pos.Reserve(m_vertex.size() + 1024);
                m_degree.Reserve(m_vertex.size() + 1024);
            }
            m_pos.set_size(m_vertex.size()+1);
            m_degree.set_size(m_vertex.size());
            num_edges = 0;

            for (auto &x: m_degree) {
                x = INVALID_ID;
            }

            m_indices = get_indices(m_vertex, _iter);
            m_mg_indices = mg.indices(m_vertex);
            assert(m_vertex.size() == m_mg_indices.size());
            assert(m_indices.size() <= _iter.size());
            m_pos[0] = 0;
            for (size_t i = 0; i < m_vertex.size(); i++){
                IdType adj_idx = m_mg_indices[i];
                m_pos[i + 1] = m_pos[i] + std::min((IdType) m_intersect.size(), m_mg->Degree(adj_idx));
            }
            est_edges = m_pos[m_vertex.size()];
            if (m_ctn.capacity() < est_edges) m_ctn.Reserve(est_edges);

            if(m_par){
                tbb::parallel_for(tbb::blocked_range<uint64_t >(0, m_indices.size()), [this](tbb::blocked_range<uint64_t> r){
                    for (uint64_t i = r.begin(); i < r.end(); ++i) {
                        IdType v_idx = m_indices[i];
                        IdType adj_idx = m_mg_indices[v_idx];
                        IdType v_id = m_vertex[v_idx];
                        IdType* start = m_ctn.begin() + m_pos[i];
                        size_t degree = m_bounded ? m_intersect.intersect(m_mg->N(adj_idx), v_id, start)
                                                  : m_intersect.intersect(m_mg->N(adj_idx), start);
                        m_degree[v_idx] = degree;
                    }
                });
            } else {
                for (IdType v_idx: m_indices) {
                    IdType adj_idx = m_mg_indices[v_idx];
                    IdType v_id = m_vertex[v_idx];
                    IdType* start = m_ctn.begin() + m_pos[v_idx];
                    size_t degree = m_bounded ? m_intersect.intersect(m_mg->N(adj_idx), v_id, start)
                                              : m_intersect.intersect(m_mg->N(adj_idx), start);
                    m_degree[v_idx] = degree;
                    num_edges += degree;
                }
            }
        }

        VertexSet<IdType> N(IdType i) override {
            if (m_degree[i] == INVALID_ID) {
                IdType v_id = m_vertex[i];
                if (m_mg != nullptr) {
                    IdType adj_idx = m_mg_indices[i];
                    IdType* start = m_ctn.begin() + m_pos[i];
                    size_t degree = m_bounded ? m_intersect.intersect(m_mg->N(adj_idx), v_id, start)
                                              : m_intersect.intersect(m_mg->N(adj_idx), start);
                    m_degree[i] = degree;
                } else {
                    IdType* start = m_ctn.begin() + m_pos[i];
                    size_t degree = m_bounded ? m_intersect.intersect(DATA_GRAPH->NBound(v_id), start)
                                              : m_intersect.intersect(DATA_GRAPH->N(v_id), start);
                    m_degree[i] = degree;
                }
            }
            return VertexSet<IdType>(Constant::EmptyID<IdType>(), m_ctn.begin() + m_pos[i], m_degree[i]);
        }

        IdType Degree(IdType i) const override {
            return m_degree[i];
        }

        ManagedContainer indices(const VertexSet<IdType>& _to_iter) const override {
            return get_indices(m_vertex, _to_iter);
        }
    };
}
#endif //MINIGRAPH_MINIGRAPH_H
