//
// Created by ubuntu on 2/2/23.
//

#ifndef MINIGRAPH_VERTEX_SET_H
#define MINIGRAPH_VERTEX_SET_H
#include <cstdint>
#include <vector>
#include <cassert>
#include <cstddef>
#include <atomic>
#include "profiler.h"
namespace minigraph {
    using IdType = uint32_t;    // support up to 4-billion number of vertexes (2^64-1 edges)
//    using IdType = uint64_t; // support up to 2^64-1 number of vertexes (2^64-1 edges)
    constexpr IdType INVALID_ID = static_cast<IdType>(-1);

    class VertexSet {
    private:
        IdType *m_data{nullptr};
        IdType m_vid{INVALID_ID};
        uint64_t m_size{0};
        bool m_pooled{false};

        class VertexSetPool {
        private:
            std::vector<IdType *> buffer_exist;
            std::vector<IdType *> buffer_avail;
        public:
            VertexSetPool() = default;

            ~VertexSetPool() {
                for (IdType *_data: buffer_exist) {
                    delete[] _data;
                }
            }

            static VertexSetPool &Get() {
                thread_local static VertexSetPool pool;
                return pool;
            };

            IdType *AllocateWorkSpace() {
                if (buffer_avail.empty()) {
                    IdType *_data = new IdType[MAX_DEGREE + 1];
                    buffer_exist.push_back(_data);
                    buffer_avail.push_back(_data);
                    TOTAL_ALLOCATED += (MAX_DEGREE + 1) * sizeof(IdType);
                }
                IdType *out = buffer_avail.back();
                buffer_avail.pop_back();
                return out;
            };

            void FreeWorkSpace(IdType *_data) {
                buffer_avail.push_back(_data);
            };
        };

    public:
        inline static uint64_t MAX_DEGREE{0};
        inline static std::atomic_uint64_t TOTAL_ALLOCATED{0};
        inline static std::shared_ptr<Profiler> profiler{nullptr};

        VertexSet() = default;

        VertexSet(IdType _vid, IdType *_data, uint64_t _size) :
                m_data{_data}, m_vid{_vid},
                m_size{_size}, m_pooled{false} {};

        // TODO add fine grained buffer management
        VertexSet(size_t capacity) : m_pooled{true} {
            m_data = static_cast<IdType *>(VertexSetPool::Get().AllocateWorkSpace());
        };

        ~VertexSet() {
            if (m_pooled) VertexSetPool::Get().FreeWorkSpace(m_data);
        };

        void swap(VertexSet &other) noexcept {
            std::swap(m_data, other.m_data);
            std::swap(m_size, other.m_size);
            std::swap(m_pooled, other.m_pooled);
            std::swap(m_vid, other.m_vid);
        };

        // reference to src / pointer copy
        VertexSet(const VertexSet &src) {
            m_data = src.m_data;
            m_size = src.m_size;
            m_vid = src.m_vid;
            m_pooled = false;
        };

        // reference to src / pointer copy
        VertexSet &operator=(const VertexSet &src) {
            if (&src != this) {
                VertexSet tmp{src};
                swap(tmp);
            }
            return *this;
        };

        VertexSet(VertexSet &&src) noexcept { swap(src); };

        VertexSet &operator=(VertexSet &&src) noexcept {
            swap(src);
            return *this;
        };

        uint64_t size() const { return m_size; };
        IdType vid() const { return m_vid; };
        IdType *begin() { return m_data; };
        IdType *end() { return m_data + m_size; };
        bool pooled() const { return m_pooled; };
        const IdType *begin() const { return m_data; };
        const IdType *end() const { return m_data + m_size; };

        inline IdType &operator[](size_t i) {
            assert(i < m_size);
            return m_data[i];
        };

        inline IdType operator[](size_t i) const {
            assert(i < m_size);
            return m_data[i];
        };

        void set_size(size_t _size) {m_size = _size;};

        inline VertexSet intersect(const VertexSet &other, IdType upper) const;
        inline VertexSet intersect(const VertexSet &other) const;
        inline size_t intersect(const VertexSet &other, IdType upper, IdType *buffer) const;
        inline size_t intersect(const VertexSet &other, IdType *buffer) const;
        inline size_t intersect_cnt(const VertexSet &other, IdType upper) const;
        inline size_t intersect_cnt(const VertexSet &other) const;
        inline VertexSet subtract(const VertexSet &other, IdType upper) const;
        inline VertexSet subtract(const VertexSet &other) const;
        inline size_t subtract_cnt(const VertexSet &other, IdType upper) const;
        inline size_t subtract_cnt(const VertexSet &other) const;
        inline VertexSet bounded(IdType upper) const;
        inline size_t bounded_cnt(IdType upper) const;
        inline VertexSet remove(IdType id) const;
        inline size_t remove_cnt(IdType id) const;
        inline VertexSet indices(const VertexSet &_vertex) const;
    };

    VertexSet VertexSet::intersect(const VertexSet &other) const {
        VertexSet out(size());
        size_t idx_l = 0, idx_r = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left == right) out[out.m_size++] = left;
        }

        profiler->add_set_comp(other.vid(), idx_l + idx_r - out.size());
        profiler->add_neb_comp(other.vid(), idx_r);

        return out;
    };

    size_t VertexSet::intersect(const VertexSet &other, IdType *buffer) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left == right) buffer[out_size++] = left;
        }

        profiler->add_mg_comp(other.vid(), idx_l + idx_r - out_size);
        profiler->add_neb_comp(other.vid(), idx_r);

        return out_size;
    };

    VertexSet VertexSet::intersect(const VertexSet &other, IdType upper) const {
        VertexSet out(size());
        size_t idx_l = 0, idx_r = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left >= upper || right >= upper) break;
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left == right) out[out.m_size++] = left;
        }

        profiler->add_set_comp(other.vid(), idx_l + idx_r - out.size());
        profiler->add_neb_comp(other.vid(), idx_r);
        return out;
    };

    size_t VertexSet::intersect(const VertexSet &other, IdType upper, IdType *buffer) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left >= upper || right >= upper) break;
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left == right) buffer[out_size++] = left;
        }

        profiler->add_mg_comp(other.vid(), idx_l + idx_r - out_size);
        profiler->add_neb_comp(other.vid(), idx_r);

        return out_size;
    };

    size_t VertexSet::intersect_cnt(const VertexSet &other) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left == right) out_size++;
        }

        profiler->add_set_comp(other.vid(), idx_l + idx_r - out_size);
        profiler->add_neb_comp(other.vid(), idx_r);

        return out_size;
    };

    size_t VertexSet::intersect_cnt(const VertexSet &other, IdType upper) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left >= upper || right >= upper) break;
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left == right) out_size++;
        }

        profiler->add_set_comp(other.vid(), idx_l + idx_r - out_size);
        profiler->add_neb_comp(other.vid(), idx_r);

        return out_size;
    };

    VertexSet VertexSet::subtract(const VertexSet &other) const {
        VertexSet out(size());
        size_t idx_l = 0, idx_r = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left < right && left != other.m_vid) out[out.m_size++] = left;
        }

        size_t subtract_out_size = out.size();

        while (idx_l < size()) {
            const IdType left = m_data[idx_l++];
            if (left != other.m_vid) out[out.m_size++] = left;
        }

        profiler->add_set_comp(other.vid(), idx_l + idx_r - subtract_out_size);
        profiler->add_neb_comp(other.vid(), idx_r);

        return out;
    };

    VertexSet VertexSet::subtract(const VertexSet &other, IdType upper) const {
        VertexSet out(size());
        size_t idx_l = 0, idx_r = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left >= upper || right >= upper) break;
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left < right && left != other.m_vid) out[out.m_size++] = left;
        }

        size_t subtract_out_size = out.size();

        while (idx_l < size()) {
            const IdType left = m_data[idx_l++];
            if (left >= upper) break;
            if (left != other.m_vid) out[out.m_size++] = left;
        }

        profiler->add_set_comp(other.vid(), idx_l + idx_r - subtract_out_size);
        profiler->add_neb_comp(other.vid(), idx_r);


        return out;
    };

    size_t VertexSet::subtract_cnt(const VertexSet &other) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left < right && left != other.m_vid) out_size++;
        }

        size_t subtract_out_size = out_size;

        while (idx_l < size()) {
            const IdType left = m_data[idx_l++];
            if (left != other.m_vid) out_size++;
        }

        profiler->add_set_comp(other.vid(), idx_l + idx_r - subtract_out_size);
        profiler->add_neb_comp(other.vid(), idx_r);

        return out_size;
    };

    size_t VertexSet::subtract_cnt(const VertexSet &other, IdType upper) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left >= upper || right >= upper) break;
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
            if (left < right && left != other.m_vid) out_size++;
        }

        size_t subtract_out_size = out_size;

        while (idx_l < size()) {
            const IdType left = m_data[idx_l++];
            if (left >= upper) break;
            if (left != other.m_vid) out_size++;
        }

        profiler->add_set_comp(other.vid(), idx_l + idx_r - subtract_out_size);
        profiler->add_neb_comp(other.vid(), idx_r);

        return out_size;
    };

    VertexSet VertexSet::bounded(IdType upper) const {
        size_t idx_l = 0;
        if (size() > 64) {
            size_t count = size();
            while (count > 0) {
                size_t it = idx_l;
                size_t step = count / 2;
                it += step;
                if (m_data[it] < upper) {
                    idx_l = ++it;
                    count -= step + 1;
                } else count = step;
            }
        } else {
            while (idx_l < size() && m_data[idx_l] < upper) idx_l++;
        }
        return VertexSet(m_vid, m_data, idx_l);
    }

    size_t VertexSet::bounded_cnt(IdType upper) const {
        size_t idx_l = 0;
        if (size() > 64) {
            size_t count = size();
            while (count > 0) {
                size_t it = idx_l;
                size_t step = count / 2;
                it += step;
                if (m_data[it] < upper) {
                    idx_l = ++it;
                    count -= step + 1;
                } else count = step;
            }
        } else {
            while (idx_l < size() && m_data[idx_l] < upper) idx_l++;
        }
        return idx_l;
    }

    VertexSet VertexSet::remove(IdType upper) const {
        size_t idx_l = 0;
        if (size() > 64) {
            size_t count = size();
            while (count > 0) {
                size_t it = idx_l;
                size_t step = count / 2;
                it += step;
                if (m_data[it] < upper) {
                    idx_l = ++it;
                    count -= step + 1;
                } else count = step;
            }
        } else {
            while (idx_l < size() && m_data[idx_l] < upper) idx_l++;
        }

        if (idx_l < m_size && m_data[idx_l] == upper) {
            VertexSet out(size());
            out.m_size = m_size - 1;
            for (size_t i = 0; i < idx_l; i++) {
                out[i] = m_data[i];
            }
            for (size_t i = idx_l; i < m_size - 1; i++) {
                out[i] = m_data[i + 1];
            }
            return out;
        };

        return VertexSet(INVALID_ID, m_data, m_size);
    }

    size_t VertexSet::remove_cnt(IdType upper) const {
        size_t idx_l = 0;
        if (size() > 64) {
            size_t count = size();
            while (count > 0) {
                size_t it = idx_l;
                size_t step = count / 2;
                it += step;
                if (m_data[it] < upper) {
                    idx_l = ++it;
                    count -= step + 1;
                } else count = step;
            }
        } else {
            while (idx_l < size() && m_data[idx_l] < upper) idx_l++;
        }

        if (idx_l < m_size && m_data[idx_l] == upper) {
            return m_size - 1;
        } else {
            return m_size;
        }
    }

    VertexSet VertexSet::indices(const VertexSet &other) const {
        VertexSet out(size());
        IdType idx_l = 0, idx_r = 0;
        while (idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if (left == right) out[out.m_size++] = idx_l;
            if (left <= right) idx_l++;
            if (right <= left) idx_r++;
        }
        return out;
    }

}
#endif //MINIGRAPH_VERTEX_SET_H
