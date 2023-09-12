//
// Created by ubuntu on 1/1/23.
//
#ifndef MINIGRAPH_PROFILER_H
#define MINIGRAPH_PROFILER_H
#include "typedef.h"
#include <atomic>
#include <cstdint>
#include <vector>
#include <string>
#include <fstream>
#include <memory>
#include <algorithm>
#include <numeric>
#include <cstring>
namespace minigraph {

    class Profiler
    {
    private:
        std::vector<std::atomic_uint64_t> VID_TO_NEB_COMP; // only consider the (mg) adj scan
        std::vector<std::atomic_uint64_t> VID_TO_SET_COMP; // adj scan + intermediate scan
        std::vector<std::atomic_uint64_t> VID_TO_MG_COMP; // mg adj scan + intermidate scan
        std::vector<std::atomic_uint64_t> VID_TO_FREQ;
        uint64_t m_psize, m_vnum;
        size_t idx(uint64_t vid) const {return m_loop * m_vnum + vid;};

    public:
        inline thread_local static uint64_t m_loop{0};

        Profiler() = default;
        Profiler(uint64_t _psize, uint64_t _v_num):
                VID_TO_SET_COMP((_psize - 1) * _v_num),
                VID_TO_MG_COMP((_psize - 1) * _v_num),
                VID_TO_NEB_COMP((_psize - 1) * _v_num),
                VID_TO_FREQ((_psize - 1) * _v_num)
        {
            m_psize = _psize;
            m_vnum = _v_num;
        };
        uint64_t total_neb_comp() {
            return std::accumulate(VID_TO_NEB_COMP.begin(), VID_TO_NEB_COMP.end(), 0llu);
        };

        uint64_t total_set_comp() {
            return std::accumulate(VID_TO_SET_COMP.begin(), VID_TO_SET_COMP.end(), 0llu);
        };

        uint64_t total_mg_comp() {
            return std::accumulate(VID_TO_MG_COMP.begin(), VID_TO_MG_COMP.end(), 0llu);
        };

        void set_cur_loop(uint64_t loop_depth){ m_loop = loop_depth;};
        void add_set_comp(uint64_t vid, uint64_t comp) { VID_TO_SET_COMP.at(idx(vid)) += comp; VID_TO_FREQ.at(idx(vid))++;};
        void add_mg_comp(uint64_t vid, uint64_t comp) { VID_TO_MG_COMP.at(idx(vid)) += comp; VID_TO_FREQ.at(idx(vid))++;};
        void add_neb_comp(uint64_t vid, uint64_t comp) {VID_TO_NEB_COMP.at(idx(vid)) += comp;};

        std::vector<uint64_t> GET_VID_TO_NEB_COMP() {
            std::vector<uint64_t > out(VID_TO_NEB_COMP.size());
            for (size_t i = 0; i < VID_TO_NEB_COMP.size(); i++) {
                out.at(i) = VID_TO_NEB_COMP.at(i);
            }
            return out;
        }
        
        std::vector<uint64_t> GET_VID_TO_SET_COMP() {
            std::vector<uint64_t > out(VID_TO_SET_COMP.size());
            for (size_t i = 0; i < VID_TO_SET_COMP.size(); i++) {
                out.at(i) = VID_TO_SET_COMP.at(i);
            }
            return out;
        };

        std::vector<uint64_t> GET_VID_TO_MG_COMP() {
            std::vector<uint64_t > out(VID_TO_MG_COMP.size());
            for (size_t i = 0; i < VID_TO_MG_COMP.size(); i++) {
                out.at(i) = VID_TO_MG_COMP.at(i);
            }
            return out;
        };

        std::vector<uint64_t> GET_VID_TO_FREQ() {
            std::vector<uint64_t > out(VID_TO_FREQ.size());
            for (size_t i = 0; i < VID_TO_FREQ.size(); i++) {
                out.at(i) = VID_TO_FREQ.at(i);
            }
            return out;
        };

        void save_loop_csv(std::string out_path) {
            constexpr char COMMA = ',';
            std::ofstream outfile;
            outfile.open(out_path + ".csv", std::ios::out);
            outfile << "ID,Loop,SetComp,GmComp,NebComp,TotalComp,ReadFreq\n";
            uint64_t max_depth = m_psize - 1;

            auto _VID_TO_SET_COMP = GET_VID_TO_SET_COMP();
            auto _VID_TO_MG_COMP = GET_VID_TO_MG_COMP();
            auto _VID_TO_NEB_COMP = GET_VID_TO_NEB_COMP();
            auto _VID_TO_FREQ = GET_VID_TO_FREQ();
            for (uint64_t loop = 0; loop < max_depth; ++loop) {
                for (uint64_t vid = 0; vid < m_vnum ; vid++) {
                    uint64_t idx = loop * m_vnum + vid;
                    uint64_t cur_set_comp =  _VID_TO_SET_COMP.at(idx) ;
                    uint64_t cur_gm_comp  =  _VID_TO_MG_COMP.at(idx) ;
                    uint64_t cur_neb_comp  =  _VID_TO_NEB_COMP.at(idx) ;
                    uint64_t cur_adj_read =  _VID_TO_FREQ.at(idx) ;
                    uint64_t total_comp   =  cur_set_comp + cur_gm_comp;
                    outfile << vid << COMMA << loop << COMMA << cur_set_comp << COMMA << cur_gm_comp << COMMA << cur_neb_comp << COMMA << total_comp << COMMA << cur_adj_read << "\n";
                }
            }
            outfile.flush();
            outfile.close();
        };

        void save_aggr_csv(std::string out_path){
            constexpr char COMMA = ',';
            std::ofstream outfile;
            outfile.open(out_path, std::ios::out);
            outfile << "ID,SetComp,GmComp,NebComp,TotalComp,ReadFreq\n";
            uint64_t max_depth = m_psize - 1;
            auto _VID_TO_SET_COMP = GET_VID_TO_SET_COMP();
            auto _VID_TO_MG_COMP = GET_VID_TO_MG_COMP();
            auto _VID_TO_NEB_COMP = GET_VID_TO_NEB_COMP();
            auto _VID_TO_FREQ = GET_VID_TO_FREQ();
            for (uint64_t vid = 0; vid < m_vnum ; vid++) {
                uint64_t cur_set_comp{0}, cur_gm_comp{0}, cur_neb_comp{0}, total_comp{0}, cur_adj_read{0};
                for (uint64_t loop = 0; loop < max_depth; ++loop) {
                    uint64_t idx = loop * m_vnum + vid;
                    cur_set_comp +=  _VID_TO_SET_COMP.at(idx) ;
                    cur_gm_comp  +=  _VID_TO_MG_COMP.at(idx) ;
                    cur_neb_comp +=  _VID_TO_NEB_COMP.at(idx);
                    cur_adj_read +=  _VID_TO_FREQ.at(idx) ;
                }
                total_comp  =  cur_set_comp + cur_gm_comp;
                outfile << vid << COMMA << cur_set_comp << COMMA << cur_gm_comp << COMMA << total_comp << COMMA << cur_adj_read << "\n";
            }
            outfile.flush();
            outfile.close();
        };

        void save_aggr_vid(std::string out_path, uint64_t * VID_TO_DEG, uint64_t * VID_TO_OFFSET, uint64_t * VID_TO_TRIANGLES) {
            constexpr char COMMA = ',';
            std::ofstream outfile;
            outfile.open(out_path + "_vid.csv", std::ios::out);
            outfile << "ID,Degree,Offset,Triangles,SetComp,GmComp,NebComp,TotalComp,ReadFreq\n";
            uint64_t max_depth = m_psize - 1;
            auto _VID_TO_SET_COMP = GET_VID_TO_SET_COMP();
            auto _VID_TO_MG_COMP = GET_VID_TO_MG_COMP();
            auto _VID_TO_NEB_COMP = GET_VID_TO_NEB_COMP();
            auto _VID_TO_FREQ = GET_VID_TO_FREQ();
            for (uint64_t vid = 0; vid < m_vnum ; vid++) {
                uint64_t cur_set_comp{0}, cur_gm_comp{0}, cur_neb_comp{0}, total_comp{0}, cur_adj_read{0};
                uint64_t deg = VID_TO_DEG[vid];
                uint64_t offset = VID_TO_OFFSET[vid];
                uint64_t triangles = VID_TO_TRIANGLES[vid];
                for (uint64_t loop = 0; loop < max_depth; ++loop) {
                    uint64_t idx = loop * m_vnum + vid;
                    cur_set_comp +=  _VID_TO_SET_COMP.at(idx) ;
                    cur_gm_comp  +=  _VID_TO_MG_COMP.at(idx) ;
                    cur_neb_comp +=  _VID_TO_NEB_COMP.at(idx);
                    cur_adj_read +=  _VID_TO_FREQ.at(idx) ;
                }
                total_comp = cur_set_comp + cur_gm_comp;
                outfile << vid          << COMMA << deg         << COMMA << offset       << COMMA << triangles    << COMMA \
                        << cur_set_comp << COMMA << cur_gm_comp << COMMA << cur_neb_comp << COMMA << total_comp   << COMMA << cur_adj_read << "\n";

            }
            outfile.flush();
            outfile.close();
        }

        void save_aggr(std::string out_path, uint64_t * VID_TO_DEG, uint64_t * VID_TO_OFFSET, uint64_t * VID_TO_TRIANGLES) {
            constexpr char COMMA = ',';
            std::ofstream degFile, offsetFile, triangleFile;
            degFile.open(out_path + "_deg.csv", std::ios::out);
            offsetFile.open(out_path + "_offset.csv", std::ios::out);
            triangleFile.open(out_path + "_triangle.csv", std::ios::out);
            degFile << "Degree,SetComp,GmComp,NebComp,TotalComp,ReadFreq\n";
            offsetFile << "Offset,SetComp,GmComp,NebComp,TotalComp,ReadFreq\n";
            triangleFile << "Triangle,SetComp,GmComp,NebComp,TotalComp,ReadFreq\n";
            uint64_t max_depth = m_psize - 1;
            auto _VID_TO_SET_COMP = GET_VID_TO_SET_COMP();
            auto _VID_TO_MG_COMP = GET_VID_TO_MG_COMP();
            auto _VID_TO_FREQ = GET_VID_TO_FREQ();
            auto _VID_TO_NEB_COMP = GET_VID_TO_NEB_COMP();

            uint64_t MAX_DEGREE = *std::max_element(VID_TO_DEG, VID_TO_DEG + m_vnum);
            uint64_t MAX_OFFSET = *std::max_element(VID_TO_OFFSET, VID_TO_OFFSET + m_vnum);
            uint64_t MAX_TRIANGLE = *std::max_element(VID_TO_TRIANGLES, VID_TO_TRIANGLES + m_vnum);

            auto _DEG_TO_SET_COMP = new uint64_t [MAX_DEGREE + 1];
            auto _DEG_TO_MG_COMP = new uint64_t [MAX_DEGREE + 1];
            auto _DEG_TO_NEB_COMP = new uint64_t [MAX_DEGREE + 1];
            auto _DEG_TO_FREQ = new uint64_t [MAX_DEGREE + 1];

            memset(_DEG_TO_FREQ, 0, sizeof(uint64_t) * (MAX_DEGREE + 1));
            memset(_DEG_TO_MG_COMP, 0, sizeof(uint64_t) * (MAX_DEGREE + 1));
            memset(_DEG_TO_SET_COMP, 0, sizeof(uint64_t) * (MAX_DEGREE + 1));
            memset(_DEG_TO_NEB_COMP, 0, sizeof(uint64_t) * (MAX_DEGREE + 1));

            auto _OFFSET_TO_SET_COMP = new uint64_t [MAX_OFFSET + 1];
            auto _OFFSET_TO_MG_COMP = new uint64_t [MAX_OFFSET + 1];
            auto _OFFSET_TO_NEB_COMP = new uint64_t [MAX_OFFSET + 1];
            auto _OFFSET_TO_FREQ = new uint64_t [MAX_OFFSET + 1];
            memset(_OFFSET_TO_SET_COMP, 0, sizeof(uint64_t) * (MAX_OFFSET + 1));
            memset(_OFFSET_TO_MG_COMP, 0, sizeof(uint64_t) * (MAX_OFFSET + 1));
            memset(_OFFSET_TO_NEB_COMP, 0, sizeof(uint64_t) * (MAX_OFFSET + 1));
            memset(_OFFSET_TO_FREQ, 0, sizeof(uint64_t) * (MAX_OFFSET + 1));

            auto _TRI_TO_SET_COMP = new uint64_t [MAX_TRIANGLE + 1];
            auto _TRI_TO_MG_COMP = new uint64_t [MAX_TRIANGLE + 1];
            auto _TRI_TO_NEB_COMP = new uint64_t [MAX_TRIANGLE + 1];
            auto _TRI_TO_FREQ = new uint64_t [MAX_TRIANGLE + 1];
            memset(_TRI_TO_SET_COMP, 0, sizeof(uint64_t) * (MAX_TRIANGLE + 1));
            memset(_TRI_TO_MG_COMP, 0, sizeof(uint64_t) * (MAX_TRIANGLE + 1));
            memset(_TRI_TO_NEB_COMP, 0, sizeof(uint64_t) * (MAX_TRIANGLE + 1));
            memset(_TRI_TO_FREQ, 0, sizeof(uint64_t) * (MAX_TRIANGLE + 1));

            for (uint64_t vid = 0; vid < m_vnum ; vid++) {
                uint64_t cur_set_comp{0}, cur_gm_comp{0}, cur_adj_read{0}, cur_neb_comp{0};
                uint64_t deg = VID_TO_DEG[vid];
                uint64_t offset = VID_TO_OFFSET[vid];
                uint64_t triangles = VID_TO_TRIANGLES[vid];
                for (uint64_t loop = 0; loop < max_depth; ++loop) {
                    uint64_t idx = loop * m_vnum + vid;
                    cur_set_comp +=  _VID_TO_SET_COMP.at(idx) ;
                    cur_gm_comp  +=  _VID_TO_MG_COMP.at(idx) ;
                    cur_neb_comp += _VID_TO_NEB_COMP.at(idx);
                    cur_adj_read +=  _VID_TO_FREQ.at(idx) ;
                }
                _DEG_TO_FREQ[deg] += cur_adj_read;
                _DEG_TO_MG_COMP[deg] += cur_gm_comp;
                _DEG_TO_SET_COMP[deg] += cur_set_comp;
                _DEG_TO_NEB_COMP[deg] += cur_neb_comp;

                _OFFSET_TO_FREQ[offset] += cur_adj_read;
                _OFFSET_TO_MG_COMP[offset] += cur_gm_comp;
                _OFFSET_TO_SET_COMP[offset] += cur_set_comp;
                _OFFSET_TO_NEB_COMP[offset] += cur_neb_comp;

                _TRI_TO_FREQ[triangles] += cur_adj_read;
                _TRI_TO_MG_COMP[triangles] += cur_gm_comp;
                _TRI_TO_SET_COMP[triangles] += cur_set_comp;
                _TRI_TO_NEB_COMP[triangles] += cur_neb_comp;
            }

            for (uint64_t i = 0; i <= MAX_DEGREE; i++) {
                uint64_t cur_set_comp = _DEG_TO_SET_COMP[i];
                uint64_t cur_mg_comp = _DEG_TO_MG_COMP[i];
                uint64_t cur_neb_comp = _DEG_TO_NEB_COMP[i];
                uint64_t cur_freq = _DEG_TO_FREQ[i];
                if (cur_set_comp == 0 && cur_freq == 0 && cur_mg_comp == 0) continue;
                degFile << i << COMMA << cur_set_comp << COMMA << cur_mg_comp << COMMA << cur_neb_comp << COMMA << cur_mg_comp + cur_set_comp << COMMA << cur_freq << "\n";
            }
            degFile.flush();
            degFile.close();

            for (uint64_t i = 0; i <= MAX_OFFSET; i++) {
                uint64_t cur_set_comp = _OFFSET_TO_SET_COMP[i];
                uint64_t cur_mg_comp = _OFFSET_TO_MG_COMP[i];
                uint64_t cur_neb_comp = _OFFSET_TO_NEB_COMP[i];
                uint64_t cur_freq = _OFFSET_TO_FREQ[i];
                if (cur_set_comp == 0 && cur_freq == 0 && cur_mg_comp == 0) continue;
                offsetFile << i << COMMA << cur_set_comp << COMMA << cur_mg_comp << COMMA << cur_neb_comp << COMMA << cur_mg_comp + cur_set_comp << COMMA << cur_freq << "\n";
            }
            offsetFile.flush();
            offsetFile.close();

            for (uint64_t i = 1; i <= MAX_TRIANGLE; i++) {
                uint64_t cur_set_comp = _TRI_TO_SET_COMP[i];
                uint64_t cur_mg_comp = _TRI_TO_MG_COMP[i];
                uint64_t cur_neb_comp = _TRI_TO_NEB_COMP[i];
                uint64_t cur_freq = _TRI_TO_FREQ[i];
                if (cur_set_comp == 0 && cur_freq == 0 && cur_mg_comp == 0) continue;
                triangleFile << i << COMMA << cur_set_comp << COMMA << cur_mg_comp << COMMA << cur_neb_comp << COMMA << cur_mg_comp + cur_set_comp << COMMA << cur_freq << "\n";
            }
            triangleFile.flush();
            triangleFile.close();

            delete[] _DEG_TO_SET_COMP;
            delete[] _DEG_TO_MG_COMP;
            delete[] _DEG_TO_FREQ;
            delete[] _DEG_TO_NEB_COMP;

            delete[] _OFFSET_TO_FREQ;
            delete[] _OFFSET_TO_MG_COMP;
            delete[] _OFFSET_TO_SET_COMP;
            delete[] _OFFSET_TO_NEB_COMP;

            delete[] _TRI_TO_FREQ;
            delete[] _TRI_TO_MG_COMP;
            delete[] _TRI_TO_SET_COMP;
            delete[] _TRI_TO_NEB_COMP;
        }

        void save_aggr_csv(std::string out_path, uint64_t * VID_TO_DEG, uint64_t * VID_TO_OFFSET, uint64_t * VID_TO_TRIANGLES){
            save_aggr_vid(out_path, VID_TO_DEG, VID_TO_OFFSET, VID_TO_TRIANGLES);
            save_aggr(out_path, VID_TO_DEG, VID_TO_OFFSET, VID_TO_TRIANGLES);
        };

    };
}
#endif //MINIGRAPH_PROFILER_H