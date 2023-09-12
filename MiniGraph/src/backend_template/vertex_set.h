//
// Created by ubuntu on 1/1/23.
//

#ifndef MINIGRAPH_VERTEX_SET_H
#define MINIGRAPH_VERTEX_SET_H
#include "common.h"
#include "stdint.h"
#include "assert.h"
#include <vector>
#include <atomic>
#include <mutex>
namespace minigraph {

    template<typename IdType>
    class VertexSet
    {
    private:
        IdType * m_data{nullptr};
        IdType m_vid{Constant::EmptyID<IdType>()};
        uint64_t m_size{0};
        bool m_pooled{false};
        class VertexSetPool;
    public:
        inline static uint64_t MAX_DEGREE{0};
        inline static std::atomic_uint64_t TOTAL_ALLOCATED{0};

        VertexSet(IdType _vid, IdType * _data, uint64_t _size):
            m_data{_data}, m_vid{_vid},
            m_size{_size}, m_pooled{false}{};

        VertexSet();
        ~VertexSet();

        void swap(VertexSet<IdType> &other) noexcept {
            std::swap(m_data, other.m_data);
            std::swap(m_size, other.m_size);
            std::swap(m_pooled, other.m_pooled);
            std::swap(m_vid, other.m_vid);
        };

        // reference to src / pointer copy
        VertexSet(const VertexSet<IdType> &src)
        {
            m_data= src.m_data;
            m_size= src.m_size;
            m_vid = src.m_vid;
            m_pooled = false;
        };

        // reference to src / pointer copy
        VertexSet &operator=(const VertexSet<IdType> &src){
            if (&src != this){
                VertexSet<IdType> tmp{src};
                swap(tmp);
            }
            return *this;
        };
        VertexSet(VertexSet<IdType> &&src) noexcept { swap(src); };
        VertexSet &operator=(VertexSet<IdType> &&src) noexcept { swap(src); return *this;};

        uint64_t size() const {return m_size;};
        IdType vid() const {return m_vid;};
        IdType* begin() {return m_data;};
        IdType* end() {return m_data + m_size;};
        bool pooled() const {return m_pooled;};
        const IdType *begin() const {return m_data;};
        const IdType *end() const {return m_data + m_size;};

        IdType &operator[](size_t i){
            assert(i < m_size);
            return m_data[i];
        };

        const IdType operator[](size_t i) const {
            assert(i < m_size);
            return m_data[i];
        };

        void set_size(size_t _size) {m_size = _size;};
        VertexSet intersect(const VertexSet& other, IdType upper) const;
        VertexSet intersect(const VertexSet& other) const;
        size_t intersect(const VertexSet& other, IdType upper, IdType * buffer) const;
        size_t intersect(const VertexSet& other, IdType * buffer) const;

        size_t intersect_cnt(const VertexSet& other, IdType upper) const;
        size_t intersect_cnt(const VertexSet& other) const;

        VertexSet subtract(const VertexSet& other, IdType upper) const;
        VertexSet subtract(const VertexSet& other) const;
        size_t subtract_cnt(const VertexSet& other, IdType upper) const;
        size_t subtract_cnt(const VertexSet& other) const;
        VertexSet bounded(IdType upper) const;
        size_t bounded_cnt(IdType upper) const;
        VertexSet remove(IdType id) const;
        size_t remove_cnt(IdType id) const;
        VertexSet indices(const VertexSet<IdType>& _vertex ) const;
    };

    // ---------- Implementation ----------

    template<typename IdType>
    class VertexSet<IdType>::VertexSetPool
    {
    private:
        std::vector<IdType *> buffer_exist;
        std::vector<IdType *> buffer_avail;
        inline static std::mutex m_mutex;

    public:
        VertexSetPool() = default;
        ~VertexSetPool() {
            for (IdType * _data : buffer_exist) {
                delete[] _data;
            }
        }
        static VertexSetPool &Get() {
            thread_local static VertexSetPool pool;
            return pool;
        };
        IdType * AllocateWorkSpace(){
            if (buffer_avail.size() == 0){
                IdType *_data = new IdType[VertexSet<IdType>::MAX_DEGREE];
                buffer_exist.push_back(_data);
                buffer_avail.push_back(_data);
                VertexSet<IdType>::TOTAL_ALLOCATED += (VertexSet<IdType>::MAX_DEGREE) * sizeof(IdType);
//                {
//                    std::lock_guard<std::mutex> lockGuard(m_mutex);
//                    LOG(INFO) << "VertexSetPool AllocateWorkSpace: " << ToReadableSize(TOTAL_ALLOCATED);
//                }
            }
            IdType * out = buffer_avail.back();
            buffer_avail.pop_back();
            return out;
        };

        void FreeWorkSpace(IdType * _data){
            buffer_avail.push_back(_data);
        };
    };

    template<typename IdType>
    VertexSet<IdType>::~VertexSet() {
        if(m_pooled) VertexSetPool::Get().FreeWorkSpace(m_data);
    }

    template<typename IdType>
    VertexSet<IdType>::VertexSet(): m_pooled{true} {
        m_data = static_cast<IdType *>(VertexSetPool::Get().AllocateWorkSpace());
    }

    template<typename IdType>
    VertexSet<IdType> VertexSet<IdType>::intersect(const VertexSet<IdType> &other) const {
        VertexSet<IdType> out;
        size_t idx_l = 0, idx_r = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left == right) out[out.m_size++] = left;
        }
        return out;
    };

    template<typename IdType>
    size_t VertexSet<IdType>::intersect(const VertexSet<IdType> &other, IdType* buffer) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left == right) buffer[out_size++] = left;
        }
        return out_size;
    };

    template<typename IdType>
    VertexSet<IdType> VertexSet<IdType>::intersect(const VertexSet<IdType> &other, IdType upper) const {
        VertexSet<IdType> out;
        size_t idx_l = 0, idx_r = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left >= upper || right >= upper) break;
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left == right) out[out.m_size++] = left;
        }
        return out;
    };

    template<typename IdType>
    size_t VertexSet<IdType>::intersect(const VertexSet<IdType> &other, IdType upper, IdType* buffer) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left >= upper || right >= upper) break;
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left == right) buffer[out_size++] = left;
        }
        return out_size;
    };

    template<typename IdType>
    size_t VertexSet<IdType>::intersect_cnt(const VertexSet<IdType> &other) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left == right) out_size++;
        }
        return out_size;
    };

    template<typename IdType>
    size_t VertexSet<IdType>::intersect_cnt(const VertexSet<IdType> &other, IdType upper) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left >= upper || right >= upper) break;
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left == right) out_size++;
        }
        return out_size;
    };

    template<typename IdType>
    VertexSet<IdType> VertexSet<IdType>::subtract(const VertexSet<IdType> &other) const {
        VertexSet<IdType> out;
        size_t idx_l = 0, idx_r = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left < right && left != other.m_vid) out[out.m_size++] = left;
        }
        while(idx_l < size()) {
            const IdType left = m_data[idx_l++];
            if (left != other.m_vid) out[out.m_size++] = left;
        }
        return out;
    };

    template<typename IdType>
    VertexSet<IdType> VertexSet<IdType>::subtract(const VertexSet<IdType> &other, IdType upper) const {
        VertexSet<IdType> out;
        size_t idx_l = 0, idx_r = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left >= upper || right >= upper) break;
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left < right && left != other.m_vid) out[out.m_size++] = left;
        }
        while(idx_l < size()) {
            const IdType left = m_data[idx_l++];
            if (left >= upper) break;
            if (left != other.m_vid) out[out.m_size++] = left;
        }
        return out;
    };

    template<typename IdType>
    size_t VertexSet<IdType>::subtract_cnt(const VertexSet<IdType> &other) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left < right && left != other.m_vid) out_size++;
        }
        while(idx_l < size()) {
            const IdType left = m_data[idx_l++];
            if (left != other.m_vid) out_size++;
        }
        return out_size;
    };

    template<typename IdType>
    size_t VertexSet<IdType>::subtract_cnt(const VertexSet<IdType> &other, IdType upper) const {
        size_t idx_l = 0, idx_r = 0, out_size = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left >= upper || right >= upper) break;
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
            if(left < right && left != other.m_vid) out_size++;
        }
        while(idx_l < size()) {
            const IdType left = m_data[idx_l++];
            if (left >= upper) break;
            if (left != other.m_vid) out_size++;
        }
        return out_size;
    };

    template<typename IdType>
    VertexSet<IdType> VertexSet<IdType>::bounded(IdType upper) const {
        size_t idx_l = 0;
        if (size() > 64){
            size_t count = size();
            while(count > 0){
                size_t it = idx_l;
                size_t step = count / 2;
                it += step;
                if (m_data[it] < upper) {
                    idx_l = ++it;
                    count -= step + 1;
                }
                else count = step;
            }
        } else {
            while(idx_l < size() && m_data[idx_l] < upper) idx_l++;
        }
        return VertexSet<IdType>(m_vid, m_data, idx_l);
    }

    template<typename IdType>
    size_t VertexSet<IdType>::bounded_cnt(IdType upper) const {
        size_t idx_l = 0;
        if (size() > 64){
            size_t count = size();
            while(count > 0){
                size_t it = idx_l;
                size_t step = count / 2;
                it += step;
                if (m_data[it] < upper) {
                    idx_l = ++it;
                    count -= step + 1;
                }
                else count = step;
            }
        } else {
            while(idx_l < size() && m_data[idx_l] < upper) idx_l++;
        }
        return idx_l;
    }

    template<typename IdType>
    VertexSet<IdType> VertexSet<IdType>::remove(IdType upper) const {
        size_t idx_l = 0;
        if (size() > 64){
            size_t count = size();
            while(count > 0){
                size_t it = idx_l;
                size_t step = count / 2;
                it += step;
                if (m_data[it] < upper) {
                    idx_l = ++it;
                    count -= step + 1;
                }
                else count = step;
            }
        } else {
            while(idx_l < size() && m_data[idx_l] < upper) idx_l++;
        }

        if (idx_l < m_size && m_data[idx_l] == upper){
            VertexSet out;
            out.m_size = m_size - 1;
            for (size_t i = 0; i < idx_l; i++){
                out[i] = m_data[i];
            }
            for (size_t i = idx_l; i < m_size - 1; i++){
                out[i] = m_data[i+1];
            }
            return out;
        };

        return VertexSet<IdType>(Constant::EmptyID<IdType>(), m_data, m_size);
    }

    template<typename IdType>
    size_t VertexSet<IdType>::remove_cnt(IdType upper) const {
        size_t idx_l = 0;
        if (size() > 64){
            size_t count = size();
            while(count > 0){
                size_t it = idx_l;
                size_t step = count / 2;
                it += step;
                if (m_data[it] < upper) {
                    idx_l = ++it;
                    count -= step + 1;
                }
                else count = step;
            }
        } else {
            while(idx_l < size() && m_data[idx_l] < upper) idx_l++;
        }

        if (idx_l < m_size && m_data[idx_l] == upper){
            return m_size - 1;
        } else {
            return m_size;
        }
    }

    template<typename IdType>
    VertexSet<IdType> VertexSet<IdType>::indices(const VertexSet<IdType>& other) const {
        VertexSet<IdType> out;
        size_t idx_l = 0, idx_r = 0;
        while(idx_l < size() && idx_r < other.size()) {
            const IdType left = m_data[idx_l];
            const IdType right = other[idx_r];
            if(left == right) out[out.m_size++] = idx_l;
            if(left <= right) idx_l++;
            if(right <= left) idx_r++;
        }
        return out;
    }
}
#endif //MINIGRAPH_VERTEX_SET_H
