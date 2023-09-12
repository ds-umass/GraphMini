//
// Created by ubuntu on 1/19/23.
//

#ifndef MINIGRAPH_BACKEND_H
#define MINIGRAPH_BACKEND_H
#include "vertex_set.h"
#include "graph.h"
#include "minigraph.h"
#include "profiler.h"
#include <cmath>
#include <omp.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/tick_count.h>
namespace minigraph {
    struct cc {
        alignas(64) long long count{0};
        cc& operator +=(long long c) {
            count += c;
            return *this;
        }
    };

    struct Context
    {
        std::shared_ptr<Profiler> profiler;
        int num_threads{1};
        int iep_redundency{1};
        tbb::tick_count tick_begin{tbb::tick_count::now()};
        std::vector<cc> per_thread_result;
        std::vector<cc> per_thread_handled;
        std::vector<double> per_thread_time; // omp
        std::vector<tbb::tick_count> per_thread_tick; // tbb
        Context(int _num_threads): num_threads{_num_threads}{
            tick_begin = tbb::tick_count::now();
            per_thread_tick.resize(num_threads, tick_begin);
            per_thread_result.resize(num_threads);
            per_thread_time.resize(num_threads);
            per_thread_handled.resize(num_threads);
        };

        double tick_time(size_t i) {
            return (per_thread_tick.at(i) - tick_begin).seconds();
        }

        long long get_handled() {
            long long out = 0;
            for (cc c: per_thread_handled) {
                out += c.count;
            }
            return out;
        }
        std::vector<size_t> get_ids() {
            // sorted by time
            std::vector<std::pair<double, int>> time_to_id;
            for (int i = 0 ; i < num_threads; i++){
                time_to_id.push_back(std::make_pair<>(per_thread_time.at(i) + (per_thread_tick.at(i) - tick_begin).seconds(), i));
            }
            std::sort(time_to_id.begin(), time_to_id.end(), [](std::pair<double, int> left, std::pair<double, int> right ) {
                return left.first < right.first;});
            std::vector<size_t> out;
            for (auto x: time_to_id) {
                out.push_back(x.second);
            }
            return out;
        }
        long long get_result(){
            long long out = 0;
            for (cc c: per_thread_result) {
                out += c.count;
            }
            return out / std::max(1, iep_redundency);
        }
        double get_max_time(){
            int i = 0;
            double out = per_thread_time.at(0) + tick_time(0);
            for (double x: per_thread_time) {
                out = std::max(x + tick_time(i), out);
                i++;
            }
            return out;
        }
        double get_min_time(){
            int i = 0;
            double out = per_thread_time.at(0) + tick_time(0);
            for (double x: per_thread_time) {
                out = std::min(x + tick_time(i), out);
                i++;
            }
            return out;
        }
        double get_mean_time() {
            int i = 0;
            double out = 0;
            for (double x: per_thread_time) {
                out += x + tick_time(i);
                i++;
            }
            return out / num_threads;
        };

        double get_var_time(){
            double mean = get_mean_time();
            double var = 0.0;
            int i = 0;
            for (double x: per_thread_time) {
                var += (x + tick_time(i) - mean) * (x + tick_time(i) - mean);
                i++;
            }
            return var / num_threads;
        };
    };
}
#endif