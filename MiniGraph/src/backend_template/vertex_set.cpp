//
// Created by ubuntu on 1/1/23.
//

#include "vertex_set.h"
namespace minigraph {
//    template<typename IdType>
//    class VertexSet<IdType>::VertexSetPool
//    {
//    private:
//        std::vector<IdType *> buffer_exist;
//        std::vector<IdType *> buffer_avail;
//        inline static std::mutex m_mutex;
//    public:
//        VertexSetPool() = default;
//        ~VertexSetPool() {
//            for (IdType * _data : buffer_exist) {
//                delete[] _data;
//            }
//        }
//        static VertexSetPool &Get() {
//            thread_local static VertexSetPool pool;
//            return pool;
//        };
//        IdType * AllocateWorkSpace(){
//            if (buffer_avail.size() == 0){
//                IdType *_data = new IdType[VertexSet<IdType>::MAX_DEGREE];
//                buffer_exist.push_back(_data);
//                buffer_avail.push_back(_data);
//                std::lock_guard<std::mutex> lock(m_mutex);
//                VertexSet<IdType>::TOTAL_ALLOCATED += (VertexSet<IdType>::MAX_DEGREE) * sizeof(IdType);
//            }
//            IdType * out = buffer_avail.back();
//            buffer_avail.pop_back();
//            return out;
//        };
//
//        void FreeWorkSpace(IdType * _data){
//            buffer_avail.push_back(_data);
//        };
//    };
//
//    template<typename IdType>
//    void VertexSet<IdType>::swap(VertexSet<IdType> &other) noexcept {
//        std::swap(m_data, other.m_data);
//        std::swap(m_size, other.m_size);
//        std::swap(m_pooled, other.m_pooled);
//        std::swap(m_vid, other.m_vid);
//    }
//
//    template<typename IdType>
//    VertexSet<IdType>::~VertexSet() {
//        if(m_pooled) VertexSetPool::Get().FreeWorkSpace(m_data);
//    }
//
//    template<typename IdType>
//    VertexSet<IdType>::VertexSet(): m_pooled{true} {
//        m_data = static_cast<IdType *>(VertexSetPool::Get().AllocateWorkSpace());
//    }
//
//    template<typename IdType>
//    VertexSet<IdType> VertexSet<IdType>::intersect(const VertexSet<IdType> &other) const {
//        VertexSet<IdType> out;
//        size_t idx_l = 0, idx_r = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left == right) out[out.m_size++] = left;
//        }
//        return out;
//    };
//
//    template<typename IdType>
//    size_t VertexSet<IdType>::intersect(const VertexSet<IdType> &other, IdType* buffer) const {
//        size_t idx_l = 0, idx_r = 0, out_size = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left == right) buffer[out_size++] = left;
//        }
//        return out_size;
//    };
//
//    template<typename IdType>
//    VertexSet<IdType> VertexSet<IdType>::intersect(const VertexSet<IdType> &other, IdType upper) const {
//        VertexSet<IdType> out;
//        size_t idx_l = 0, idx_r = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left >= upper || right >= upper) break;
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left == right) out[out.m_size++] = left;
//        }
//        return out;
//    };
//
//    template<typename IdType>
//    size_t VertexSet<IdType>::intersect(const VertexSet<IdType> &other, IdType upper, IdType* buffer) const {
//        size_t idx_l = 0, idx_r = 0, out_size = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left >= upper || right >= upper) break;
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left == right) buffer[out_size++] = left;
//        }
//        return out_size;
//    };
//
//    template<typename IdType>
//    size_t VertexSet<IdType>::intersect_cnt(const VertexSet<IdType> &other) const {
//        size_t idx_l = 0, idx_r = 0, out_size = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left == right) out_size++;
//        }
//        return out_size;
//    };
//
//    template<typename IdType>
//    size_t VertexSet<IdType>::intersect_cnt(const VertexSet<IdType> &other, IdType upper) const {
//        size_t idx_l = 0, idx_r = 0, out_size = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left >= upper || right >= upper) break;
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left == right) out_size++;
//        }
//        return out_size;
//    };
//
//    template<typename IdType>
//    VertexSet<IdType> VertexSet<IdType>::subtract(const VertexSet<IdType> &other) const {
//        VertexSet<IdType> out;
//        size_t idx_l = 0, idx_r = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left < right && left != other.m_vid) out[out.m_size++] = left;
//        }
//        while(idx_l < size()) {
//            const IdType left = m_data[idx_l++];
//            if (left != other.m_vid) out[out.m_size++] = left;
//        }
//        return out;
//    };
//
//    template<typename IdType>
//    VertexSet<IdType> VertexSet<IdType>::subtract(const VertexSet<IdType> &other, IdType upper) const {
//        VertexSet<IdType> out;
//        size_t idx_l = 0, idx_r = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left >= upper || right >= upper) break;
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left < right && left != other.m_vid) out[out.m_size++] = left;
//        }
//        while(idx_l < size()) {
//            const IdType left = m_data[idx_l++];
//            if (left >= upper) break;
//            if (left != other.m_vid) out[out.m_size++] = left;
//        }
//        return out;
//    };
//
//    template<typename IdType>
//    size_t VertexSet<IdType>::subtract_cnt(const VertexSet<IdType> &other) const {
//        size_t idx_l = 0, idx_r = 0, out_size = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left < right && left != other.m_vid) out_size++;
//        }
//        while(idx_l < size()) {
//            const IdType left = m_data[idx_l++];
//            if (left != other.m_vid) out_size++;
//        }
//        return out_size;
//    };
//
//    template<typename IdType>
//    size_t VertexSet<IdType>::subtract_cnt(const VertexSet<IdType> &other, IdType upper) const {
//        size_t idx_l = 0, idx_r = 0, out_size = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left >= upper || right >= upper) break;
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left < right && left != other.m_vid) out_size++;
//        }
//        while(idx_l < size()) {
//            const IdType left = m_data[idx_l++];
//            if (left >= upper) break;
//            if (left != other.m_vid) out_size++;
//        }
//        return out_size;
//    };
//
//    template<typename IdType>
//    VertexSet<IdType> VertexSet<IdType>::bounded(IdType upper) const {
//        size_t idx_l = 0;
//        if (size() > 64){
//            size_t count = size();
//            while(count > 0){
//                size_t it = idx_l;
//                size_t step = count / 2;
//                it += step;
//                if (m_data[it] < upper) {
//                    idx_l = ++it;
//                    count -= step + 1;
//                }
//                else count = step;
//            }
//        } else {
//            while(idx_l < size() && m_data[idx_l] < upper) idx_l++;
//        }
//        return VertexSet<IdType>(m_vid, m_data, idx_l);
//    }
//
//    template<typename IdType>
//    size_t VertexSet<IdType>::bounded_cnt(IdType upper) const {
//        size_t idx_l = 0;
//        if (size() > 64){
//            size_t count = size();
//            while(count > 0){
//                size_t it = idx_l;
//                size_t step = count / 2;
//                it += step;
//                if (m_data[it] < upper) {
//                    idx_l = ++it;
//                    count -= step + 1;
//                }
//                else count = step;
//            }
//        } else {
//            while(idx_l < size() && m_data[idx_l] < upper) idx_l++;
//        }
//        return idx_l;
//    }
//
//    template<typename IdType>
//    VertexSet<IdType> VertexSet<IdType>::remove(IdType upper) const {
//        size_t idx_l = 0;
//        if (size() > 64){
//            size_t count = size();
//            while(count > 0){
//                size_t it = idx_l;
//                size_t step = count / 2;
//                it += step;
//                if (m_data[it] < upper) {
//                    idx_l = ++it;
//                    count -= step + 1;
//                }
//                else count = step;
//            }
//        } else {
//            while(idx_l < size() && m_data[idx_l] < upper) idx_l++;
//        }
//
//        if (idx_l < m_size && m_data[idx_l] == upper){
//            VertexSet out;
//            out.m_size = m_size - 1;
//            for (size_t i = 0; i < idx_l; i++){
//                out[i] = m_data[i];
//            }
//            for (size_t i = idx_l; i < m_size - 1; i++){
//                out[i] = m_data[i+1];
//            }
//            return out;
//        };
//
//        return VertexSet<IdType>(Constant::EmptyID<IdType>(), m_data, m_size);
//    }
//
//    template<typename IdType>
//    size_t VertexSet<IdType>::remove_cnt(IdType upper) const {
//        size_t idx_l = 0;
//        if (size() > 64){
//            size_t count = size();
//            while(count > 0){
//                size_t it = idx_l;
//                size_t step = count / 2;
//                it += step;
//                if (m_data[it] < upper) {
//                    idx_l = ++it;
//                    count -= step + 1;
//                }
//                else count = step;
//            }
//        } else {
//            while(idx_l < size() && m_data[idx_l] < upper) idx_l++;
//        }
//
//        if (idx_l < m_size && m_data[idx_l] == upper){
//            return m_size - 1;
//        } else {
//            return m_size;
//        }
//    }
//
//    template<typename IdType>
//    VertexSet<IdType> VertexSet<IdType>::indices(const VertexSet<IdType>& other) const {
//        VertexSet<IdType> out;
//        size_t idx_l = 0, idx_r = 0;
//        while(idx_l < size() && idx_r < other.size()) {
//            const IdType left = m_data[idx_l];
//            const IdType right = other[idx_r];
//            if(left <= right) idx_l++;
//            if(right <= left) idx_r++;
//            if(left == right) out[out.m_size++] = left;
//        }
//        return out;
//    }
//
    template class VertexSet<uint64_t >;
    template class VertexSet<uint32_t >;
}

