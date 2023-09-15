#include "../include/schedule.h"
#include "../include/graph.h"
#include "../include/dataloader.h"
#include <cstdio>
#include <cstring>
#include <assert.h>
#include <algorithm>

void print(std::vector<int> vec) {
    for (auto i: vec) {
        printf("%d ", i);
    }
}

class RankClass {
private:
    std::vector<int> m_order{};
    std::vector<int> m_in_degrees{};
    std::vector<int> m_weighted_in_degrees{};
    std::vector<int> m_out_degrees{};
    std::vector<int> m_rank{};
public:
    RankClass() = default;

    RankClass(std::vector<int> order,
              std::vector<int> in_degrees,
              std::vector<int> weighted_in_degrees,
              std::vector<int> out_degrees) :
            m_order{order},
            m_in_degrees{in_degrees},
            m_weighted_in_degrees{weighted_in_degrees},
            m_out_degrees{out_degrees} {
        m_rank.resize(m_order.size());
        for (int i = 0; i < m_order.size(); i++) m_rank[m_order[i]] = i;
    };

    std::vector<int> get_order() const { return m_order; };

    void print() {
        printf("rank: ");
        for (auto r: m_order) {
            printf("%d ", r);
        }

        printf("\nInDegree:\n");
        for (auto d: m_in_degrees) {
            printf("%d ", d);
        }
        printf("\nWeightedInDegree:\n");
        for (auto d: m_weighted_in_degrees) {
            printf("%d ", d);
        }

        printf("\nOutDegree:\n");
        for (auto d: m_out_degrees) {
            printf("%d ", d);
        }
        printf("\n");
    }

    void print(std::vector<int> &adj_mat, int size) {
        std::vector<int> reordered_adj_mat;
        reordered_adj_mat.resize(size * size);

        printf("%d ", size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                reordered_adj_mat[INDEX(m_rank[i], m_rank[j], size)] = adj_mat[INDEX(i, j, size)];
            }
        }

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                printf("%d", reordered_adj_mat[INDEX(i, j, size)]);
            }
            printf(" ");
        }
        printf("\n");
    }

    // for sorting
    bool operator<(RankClass const &other) const {
        for (int i = 0; i < m_order.size(); i++) {
            if (m_in_degrees.at(i) > other.m_in_degrees.at(i)) return true;
            if (m_in_degrees.at(i) < other.m_in_degrees.at(i)) return false;
        }


        for (int i = 0; i < m_order.size(); i++) {
            if (m_weighted_in_degrees.at(i) > other.m_weighted_in_degrees.at(i)) return true;
            if (m_weighted_in_degrees.at(i) < other.m_weighted_in_degrees.at(i)) return false;
        }

        for (int i = 0; i < m_order.size(); i++) {
            if (m_out_degrees.at(i) > other.m_out_degrees.at(i)) return true;
            if (m_out_degrees.at(i) < other.m_out_degrees.at(i)) return false;
        }

        for (int i = 0; i < m_order.size(); i++) {
            if (m_order.at(i) > other.m_order.at(i)) return true;
            if (m_order.at(i) < other.m_order.at(i)) return false;
        }
        return false;
    }

};

class PatternGraph {
private:
    std::vector<int> m_adj_mat;
    int m_size;

    void get_full_permutation(std::vector<std::vector<int> > &vec, std::vector<bool> &use, std::vector<int> tmp_vec,
                              int depth) const {
        if (depth == m_size) {
            vec.push_back(tmp_vec);
            return;
        }
        for (int i = 0; i < m_size; ++i)
            if (use[i] == false) {
                use[i] = true;
                tmp_vec.push_back(i);
                get_full_permutation(vec, use, tmp_vec, depth + 1);
                tmp_vec.pop_back();
                use[i] = false;
            }
    }

    void
    get_order_greedy(std::vector<std::vector<int> > &vec, std::vector<bool> &use, std::vector<int> tmp_vec, int depth) {
        if (depth == m_size) {
            assert(tmp_vec.size() == depth);
            vec.push_back(tmp_vec);
            return;
        }
        int max_in_degree = 0;
        for (int i = 0; i < m_size; ++i) {
            if (use[i]) continue;
            int in_degree = InDegree(i, use);
            if (max_in_degree < in_degree) {
                max_in_degree = in_degree;
            }
        }

        std::vector<int> candidates;
        for (int i = 0; i < m_size; ++i) {
            if (use[i]) continue;
            if (max_in_degree == InDegree(i, use)) {
                candidates.push_back(i);
            }
        }

        int max_out_degree = 0;
        for (int i: candidates) {
            int out_degree = OutDegree(i, use);
            if (max_out_degree < out_degree) {
                max_out_degree = out_degree;
            }
        }
        std::vector<int> new_candidates;
        for (int i: candidates) {
            if (max_out_degree == OutDegree(i, use)) {
                new_candidates.push_back(i);
            }
        }

        for (int i: new_candidates) {
            assert(use[i] == false);
            use[i] = true;
            tmp_vec.push_back(i);
            get_order_greedy(vec, use, tmp_vec, depth + 1);
            tmp_vec.pop_back();
            use[i] = false;
        }

    }

    int OutDegree(int v_idx, const std::vector<bool> &use) {
        int res = 0;
        for (int j = 0; j < m_size; ++j) {
            if (use[j] == false && j != v_idx) {
                res += m_adj_mat[INDEX(v_idx, j, m_size)];
            }
        }
        return res;
    }

    int InDegree(int v_idx, const std::vector<bool> &use) {
        int res = 0;
        for (int j = 0; j < m_size; ++j) {
            if (use[j] == true && j != v_idx) res += m_adj_mat[INDEX(j, v_idx, m_size)];
        }
        return res;
    }

    int pow2(int pow) {
        int result = 1;
        for (int i = 1; i <= pow; i++) result *= 2;
        return result;
    }

    int WeightedInDegree(int v_idx, const std::vector<bool> &use, const std::vector<int> &rank) {
        int res = 0;
        for (int j = 0; j < m_size; ++j) {
            if (use[j] == true) {
                auto j_rank = std::find(rank.begin(), rank.end(), j) - rank.begin();
                res += m_adj_mat[INDEX(j, v_idx, m_size)] * pow2(m_size - j_rank);
            }
        }
        return res;
    };

    std::vector<int> InDegrees(const std::vector<int> &rank) {
        std::vector<int> result;
        std::vector<bool> use;
        use.resize(m_size);
        for (auto r: rank) {
            use[r] = true;
            result.push_back(InDegree(r, use));
        }
        return result;
    }

    std::vector<int> WeightedInDegrees(const std::vector<int> &rank) {
        std::vector<int> result;
        std::vector<bool> use;
        use.resize(m_size);
        for (auto r: rank) {
            use[r] = true;
            result.push_back(WeightedInDegree(r, use, rank));
        }
        return result;
    }

    std::vector<int> OutDegrees(const std::vector<int> &rank) {
        std::vector<int> result;
        std::vector<bool> use;
        use.resize(m_size);
        for (auto r: rank) {
            use[r] = true;
            result.push_back(OutDegree(r, use));
        }
        return result;
    }

public:
    PatternGraph(const int *adj_mat, int size) {
        m_size = size;
        m_adj_mat.clear();
        m_adj_mat.resize(m_size * m_size);
        for (int i = 0; i < m_size * m_size; i++) {
            m_adj_mat[i] = adj_mat[i];
        }

        // remove self loop
        for (int i = 0; i < m_size; i++) {
            m_adj_mat[INDEX(i, i, m_size)] = 0;
        }

        // check is undirected
        for (int i = 0; i < m_size; i++) {
            for (int j = 0; j < m_size; j++) {
                assert(m_adj_mat[INDEX(i, j, m_size)] == m_adj_mat[INDEX(j, i, m_size)]);
            }
        }
    }

    PatternGraph(const char *adj_mat, int size) {
        m_size = size;
        m_adj_mat.clear();
        m_adj_mat.resize(m_size * m_size);
        for (int i = 0; i < m_size * m_size; i++) {
            if (adj_mat[i] == '0') {
                m_adj_mat[i] = 0;
            } else {
                m_adj_mat[i] = 1;
            }
        }

        // remove self loop
        for (int i = 0; i < m_size; i++) {
            m_adj_mat[INDEX(i, i, m_size)] = 0;
        }

        // check is undirected
        for (int i = 0; i < m_size; i++) {
            for (int j = 0; j < m_size; j++) {
                assert(m_adj_mat[INDEX(i, j, m_size)] == m_adj_mat[INDEX(j, i, m_size)]);
            }
        }
    }

    std::vector<int> const get_adj_mat() { return m_adj_mat; };

    std::vector<int> get_order() {
        std::vector<std::vector<int>> ranks;
        std::vector<bool> use;
        use.resize(m_size);
        std::vector<int> rank;
//        get_order_greedy(ranks, use, rank, 0);
        get_full_permutation(ranks, use, rank, 0);

        std::vector<RankClass> rank_class_vec;
        for (int i = 0; i < (int) ranks.size(); i++) {
            auto cur_rank = ranks[i];
            auto in_degrees = InDegrees(cur_rank);
            auto weighted_in_degrees = WeightedInDegrees(cur_rank);
            auto out_degrees = OutDegrees(cur_rank);
            rank_class_vec.push_back(RankClass(cur_rank, in_degrees, weighted_in_degrees, out_degrees));
        }
        std::sort(rank_class_vec.begin(), rank_class_vec.end());
//        int idx = 0;
//        for (auto rank_class: rank_class_vec) {
//            printf("---%d---\n", idx++);
//            rank_class.print();
//            rank_class.print(m_adj_mat, m_size);
//            printf("\n");
//        }
        return rank_class_vec[0].get_order();
    };

    // return pat
    std::vector<int> get_reordered_adj_mat() {
        auto order = get_order();
        std::vector<int> rank;
        rank.resize(m_size);
        for (int i = 0; i < m_size; ++i) rank[order[i]] = i;

        std::vector<int> reordered_adj_mat;
        reordered_adj_mat.resize(m_size * m_size);
        for (int i = 0; i < m_size; i++) {
            for (int j = 0; j < m_size; j++) {
                reordered_adj_mat[INDEX(rank[i], rank[j], m_size)] = m_adj_mat[INDEX(i, j, m_size)];
            }
        }
        return reordered_adj_mat;
    };

};

void Schedule::get_greedy_schedule(const char *_adj_mat, int _size) {
    PatternGraph myPattern(_adj_mat, _size);
    auto reordered_adj_mat = myPattern.get_reordered_adj_mat();
    size = _size;
    adj_mat = new int[size * size];
    for (int i = 0; i < size * size; i++) {
        adj_mat[i] = reordered_adj_mat[i];
    }
    std::vector<std::vector<std::pair<int, int> > > restricts_vector;
    restricts_vector.clear();
    restricts_generate(adj_mat, restricts_vector);

    if (restricts_vector.size() > 0) {
        restrict_pair = restricts_vector[0];
    }

    print_schedule();
}

std::string Schedule::get_adj_mat_str() {
    std::string adj_mat_str;
    adj_mat_str.resize(size * size, '0');
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (adj_mat[INDEX(i, j, size)] == 0) {
                adj_mat_str.at(INDEX(i, j, size)) = '0';

            } else {
                adj_mat_str.at(INDEX(i, j, size)) = '1';

            }

        }
    }
    return adj_mat_str;
}

std::vector<std::vector<int>> remove_automorphisms(const char * adj_mat, int size, std::vector<std::vector<int>> orders){
    std::vector<std::string> all_mat;
    for (auto &order: orders) {
        int rank[size];
        for (int i = 0; i < size; ++i) rank[order[i]] = i;
        std::string loc_adj_mat_str;
        loc_adj_mat_str.resize(size * size);
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                loc_adj_mat_str.at(rank[i] * size + rank[j]) = adj_mat[i * size + j];
        all_mat.push_back(loc_adj_mat_str);
    }
    std::vector<std::string> org_all_mat = all_mat;
    // remove identical pattern matrix caused by automorphisms
    std::sort(all_mat.begin(), all_mat.end(), std::greater<std::string>());
    all_mat.erase(std::unique(all_mat.begin(), all_mat.end()), all_mat.end());

    std::vector<std::vector<int>> res;
    for (auto & unique_mat: all_mat){
        auto itr = std::find(org_all_mat.begin(), org_all_mat.end(), unique_mat);
        assert(itr != org_all_mat.end() && "org_all_mat does not contain unique_mat");
        auto idx = std::distance(org_all_mat.begin(), itr);
        res.push_back(orders.at(idx));
    }
    return res;
};

// GraphPi's algorithm for generating schedule in for edge-induced pattern
// This implementation additionally removes automorphisms before computing the scores for each iteration
void Schedule::get_schedule(const char *_adj_mat, int _size, uint64_t v_cnt, uint64_t e_cnt, uint64_t tri_cnt, PerfModelType model_type) {
    size = _size;
    adj_mat = new int[size * size];
    int org_adj_mat[size * size];
    for (int i = 0; i < size * size; i++) {
        if (_adj_mat[i] == '1') {
            adj_mat[i] = 1;
            org_adj_mat[i] = 1;
        } else if (_adj_mat[i] == '0') {
            adj_mat[i] = 0;
            org_adj_mat[i] = 0;
        }
    }
    std::vector<std::pair<int, int> > best_pairs;
    best_pairs.clear();
    unsigned int pow = 1;
    for (int i = 2; i <= size; ++i) pow *= i;

    std::vector<std::vector<int> > candidate_permutations;
    candidate_permutations.clear();

    bool use[size];
    for (int i = 0; i < size; ++i) use[i] = false;
    std::vector<int> tmp_vec;

    get_full_permutation(candidate_permutations, use, tmp_vec, 0);
    assert(candidate_permutations.size() == pow);
    // remove automorphism
    candidate_permutations = remove_automorphisms(_adj_mat, size, candidate_permutations);
    remove_invalid_permutation(candidate_permutations);
    int max_val = 0;
    for (const auto &vec: candidate_permutations) {
        max_val = std::max(max_val, get_vec_optimize_num(vec));
    }
    std::vector<std::vector<int> > tmp_candidate;
    tmp_candidate.clear();
    for (const auto &vec: candidate_permutations) {
        if (get_vec_optimize_num(vec) == max_val) {
            tmp_candidate.push_back(vec);
        }
    }

    candidate_permutations = tmp_candidate;

    int *best_order = new int[size];
    double min_val = 0.0;
    bool have_best = false;

    // brute force search
    for (const auto &vec: candidate_permutations) {
        int rank[size];
        for (int i = 0; i < size; ++i) rank[vec[i]] = i;

        int cur_adj_mat[size * size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector<std::vector<std::pair<int, int> > > restricts_vector;

        restricts_vector.clear();

        restricts_generate(cur_adj_mat, restricts_vector);

        if (restricts_vector.size() == 0) {
            std::vector<std::pair<int, int> > Empty;
            Empty.clear();

            double val;
            if (model_type == PerfModelType::graphpi) val = our_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt, tri_cnt);
            else val = GraphZero_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
            if (have_best == false || val < min_val) {
                have_best = true;
                min_val = val;
                for (int i = 0; i < size; ++i) best_order[i] = vec[i];
                best_pairs = Empty;
            }
        }

        for (const auto &pairs: restricts_vector) {
            double val;
            if (model_type == PerfModelType::graphpi) val = our_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt, tri_cnt);
            else val = GraphZero_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
            if (have_best == false || val < min_val) {
//                printf("cost=%.1f\n", val);
                have_best = true;
                min_val = val;
                for (int i = 0; i < size; ++i) best_order[i] = vec[i];
                best_pairs = pairs;
            }
        }

    }

    int rank[size];
    for (int i = 0; i < size; ++i) rank[best_order[i]] = i;

//    const int* pattern_adj_mat = pattern.get_adj_mat_ptr();
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            adj_mat[INDEX(rank[i], rank[j], size)] = org_adj_mat[INDEX(i, j, size)];
    delete[] best_order;

    restrict_pair = best_pairs;

    std::vector<int> I;
    I.clear();
    for (int i = 0; i < size; ++i) I.push_back(i);
    in_exclusion_optimize_num = get_vec_optimize_num(I);
    if (in_exclusion_optimize_num <= 1) {
        // printf("Can not use in_exclusion_optimize with this schedule\n");
        in_exclusion_optimize_num = 0;
    } else {
        in_exclusion_optimize_num--; // Fix a bug in the later redundancy computation
        // printf("use in_exclusion_optimize with size %d\n", in_exclusion_optimize_num);
        init_in_exclusion_optimize();
    }
    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;
    bool is_pattern_valid = true;
    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i) {
        bool valid = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)]) {
                valid = true;
                break;
            }
        if (valid == false) {
            //Invalid Schedule
            is_pattern_valid = false;
            return;
        }
    }
    assert(is_pattern_valid && 'Invalid pattern');
    build_loop_invariant();
    add_restrict(best_pairs);
    int tmp = get_in_exclusion_optimize_num();
    if (tmp <= 0) {
        in_exclusion_optimize_redundancy = 1;
    } else {
        Graph *complete;
        DataLoader *D = new DataLoader();
        assert(D->load_complete(complete, get_size()));
        delete D;
        in_exclusion_optimize_redundancy = 1;
        long long ans = complete->pattern_matching(*this, 1);
        set_in_exclusion_optimize_num(0);
        long long true_ans = complete->pattern_matching(*this, 1);
        set_in_exclusion_optimize_num(tmp + 1); // change the number back to the correct one
        init_in_exclusion_optimize();
        delete complete;
        in_exclusion_optimize_redundancy = ans / true_ans;
    }
    print_schedule();
};

Schedule::Schedule(const Pattern &pattern, bool &is_pattern_valid, int performance_modeling_type, int restricts_type,
                   bool use_in_exclusion_optimize, long long int v_cnt, long long int e_cnt, long long tri_cnt) {
    if (performance_modeling_type != 0 && tri_cnt == -1) {
        printf("Fatal: Can not use performance modeling if not have triangle number of this dataset.\n");
        fflush(stdout);
        assert(0);
    }

    is_pattern_valid = true;
    size = pattern.get_size();
    adj_mat = new int[size * size];

    // not use performance_modeling, simply copy the pat from pattern
    memcpy(adj_mat, pattern.get_adj_mat_ptr(), size * size * sizeof(int));

    std::vector<std::pair<int, int> > best_pairs;
    best_pairs.clear();
    //Initialize pat
    //If we use performance_modeling, we may change the order of vertex,
    //the best order produced by performance_modeling(...) is saved in best_order[]
    //Finally, we use best_order[] to relocate pat
    if (performance_modeling_type != 0) {
        unsigned int pow = 1;
        for (int i = 2; i <= size; ++i) pow *= i;

        std::vector<std::vector<int> > candidate_permutations;
        candidate_permutations.clear();

        bool use[size];
        for (int i = 0; i < size; ++i) use[i] = false;
        std::vector<int> tmp_vec;
        get_full_permutation(candidate_permutations, use, tmp_vec, 0);
        assert(candidate_permutations.size() == pow);

        remove_invalid_permutation(candidate_permutations);

        if (performance_modeling_type == 1) {
            //reduce candidates
            int max_val = 0;
            for (const auto &vec: candidate_permutations) {
                max_val = std::max(max_val, get_vec_optimize_num(vec));
            }
            std::vector<std::vector<int> > tmp;
            tmp.clear();
            for (const auto &vec: candidate_permutations)
                if (get_vec_optimize_num(vec) == max_val) {
                    tmp.push_back(vec);
                }
            candidate_permutations = tmp;
        }

        int *best_order = new int[size];
        double min_val;
        bool have_best = false;


        for (const auto &vec: candidate_permutations) {
            int rank[size];
            for (int i = 0; i < size; ++i) rank[vec[i]] = i;

            int *cur_adj_mat;
            cur_adj_mat = new int[size * size];
            for (int i = 0; i < size; ++i)
                for (int j = 0; j < size; ++j)
                    cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

            std::vector<std::vector<std::pair<int, int> > > restricts_vector;

            restricts_vector.clear();

            if (restricts_type == 1) {
                restricts_generate(cur_adj_mat, restricts_vector);
            } else {
                Schedule schedule(cur_adj_mat, size);

                std::vector<std::pair<int, int> > pairs;
                schedule.GraphZero_aggressive_optimize(pairs);

                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }

            if (restricts_vector.size() == 0) {
                std::vector<std::pair<int, int> > Empty;
                Empty.clear();

                double val;
                if (performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt, tri_cnt);
                } else {
                    if (performance_modeling_type == 2) {
                        val = GraphZero_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    } else {
                        val = Naive_estimate_schedule_restrict(vec, Empty, v_cnt, e_cnt);
                    }
                }

                if (have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for (int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = Empty;
                }
            }


            for (const auto &pairs: restricts_vector) {
                double val;
                if (performance_modeling_type == 1) {
                    val = our_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt, tri_cnt);
                } else {
                    if (performance_modeling_type == 2) {
                        val = GraphZero_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                    } else {
                        val = Naive_estimate_schedule_restrict(vec, pairs, v_cnt, e_cnt);
                    }
                }

                if (have_best == false || val < min_val) {
                    have_best = true;
                    min_val = val;
                    for (int i = 0; i < size; ++i) best_order[i] = vec[i];
                    best_pairs = pairs;
                }
            }

        }

        int rank[size];
        for (int i = 0; i < size; ++i) rank[best_order[i]] = i;

        const int *pattern_adj_mat = pattern.get_adj_mat_ptr();
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                adj_mat[INDEX(rank[i], rank[j], size)] = pattern_adj_mat[INDEX(i, j, size)];
        delete[] best_order;
    } else {
        std::vector<int> I;
        I.clear();
        for (int i = 0; i < size; ++i) I.push_back(i);

        std::vector<std::vector<std::pair<int, int> > > restricts_vector;
        restricts_vector.clear();

        if (restricts_type != 0) {

            if (restricts_type == 1) {
                restricts_generate(adj_mat, restricts_vector);
            } else {
                std::vector<std::pair<int, int> > pairs;
                GraphZero_aggressive_optimize(pairs);

                restricts_vector.clear();
                restricts_vector.push_back(pairs);
            }
        }

        bool have_best = false;
        double min_val;

        for (const auto &pairs: restricts_vector) {
            double val;
            if (restricts_type == 1) {
                val = our_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt, tri_cnt);
            } else {
                val = GraphZero_estimate_schedule_restrict(I, pairs, v_cnt, e_cnt);
            }
            if (have_best == false || val < min_val) {
                have_best = true;
                min_val = val;
                best_pairs = pairs;
            }
        }

    }

    if (use_in_exclusion_optimize) {
        std::vector<int> I;
        I.clear();
        for (int i = 0; i < size; ++i) I.push_back(i);
        in_exclusion_optimize_num = get_vec_optimize_num(I);
        if (in_exclusion_optimize_num <= 1) {
            printf("Can not use in_exclusion_optimize with this schedule\n");
            in_exclusion_optimize_num = 0;
        } else {
            printf("use in_exclusion_optimize with size %d\n", in_exclusion_optimize_num);
            init_in_exclusion_optimize();
        }
    } else {
        in_exclusion_optimize_num = 0;
    }

    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;

    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i) {
        bool valid = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)]) {
                valid = true;
                break;
            }
        if (valid == false) {
            //Invalid Schedule
            is_pattern_valid = false;
            return;
        }
    }

    build_loop_invariant();
    if (restricts_type != 0) add_restrict(best_pairs);

    set_in_exclusion_optimize_redundancy();
}

Schedule::Schedule(const int *_adj_mat, int _size) {
    size = _size;
    adj_mat = new int[size * size];

    memcpy(adj_mat, _adj_mat, size * size * sizeof(int));

    // The I-th loop consists of at most the intersection of i-1 VertexSet.
    // So the max number of prefix = 0 + 1 + ... + size-1 = size * (size-1) / 2
    int max_prefix_num = size * (size - 1) / 2;
    father_prefix_id = new int[max_prefix_num];
    last = new int[size];
    next = new int[max_prefix_num];
    loop_set_prefix_id = new int[size];
    prefix = new Prefix[max_prefix_num];
    restrict_last = new int[size];
    restrict_next = new int[max_prefix_num];
    restrict_index = new int[max_prefix_num];
    memset(father_prefix_id, -1, max_prefix_num * sizeof(int));
    memset(last, -1, size * sizeof(int));
    memset(next, -1, max_prefix_num * sizeof(int));
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));

    total_prefix_num = 0;
    total_restrict_num = 0;
    in_exclusion_optimize_num = 0;

    // The I-th vertex must connect with at least one vertex from 0 to i-1.
    for (int i = 1; i < size; ++i) {
        bool valid = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)]) {
                valid = true;
                break;
            }
        if (valid == false) {
            printf("invalid schedule!\n");
            assert(0);
        }
    }

    build_loop_invariant();

    set_in_exclusion_optimize_redundancy();
}

Schedule::~Schedule() {
    delete[] adj_mat;
    delete[] father_prefix_id;
    delete[] last;
    delete[] next;
    delete[] loop_set_prefix_id;
    delete[] prefix;
    delete[] restrict_last;
    delete[] restrict_next;
    delete[] restrict_index;
}

int Schedule::get_max_degree() const {
    int mx = 0;
    for (int i = 0; i < size; ++i) {
        int cnt = 0;
        for (int j = 0; j < size; ++j)
            cnt += adj_mat[INDEX(i, j, size)];
        if (cnt > mx) mx = cnt;
    }
    return mx;
}

void Schedule::build_loop_invariant() {
    int *tmp_data = new int[size];
    loop_set_prefix_id[0] = -1;
    for (int i = 1; i < size; ++i) {
        int data_size = 0;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(i, j, size)])
                tmp_data[data_size++] = j;
        loop_set_prefix_id[i] = find_father_prefix(data_size, tmp_data);
    }
    assert(total_prefix_num <= size * (size - 1) / 2);
    delete[] tmp_data;
}

int Schedule::find_father_prefix(int data_size, const int *data) {
    if (data_size == 0)
        return -1;
    int num = data[data_size - 1];
    for (int prefix_id = last[num]; prefix_id != -1; prefix_id = next[prefix_id])
        if (prefix[prefix_id].equal(data_size, data))
            return prefix_id;

    // not found, create new prefix and find its father prefix id recursively
    int father = find_father_prefix(data_size - 1, data);
    father_prefix_id[total_prefix_num] = father;
    next[total_prefix_num] = last[num];
    last[num] = total_prefix_num;
    prefix[total_prefix_num].init(data_size, data);
    ++total_prefix_num;
    return total_prefix_num - 1;
}

void Schedule::add_restrict(const std::vector<std::pair<int, int> > &restricts) {
    restrict_pair = restricts;
    for (unsigned int i = 0; i < restrict_pair.size();) {
        bool tag = true;
        for (unsigned int j = 0; j < restrict_pair.size(); ++j) {
            if (i != j && restrict_pair[j].first == restrict_pair[i].first)
                for (unsigned int k = 0; k < restrict_pair.size(); ++k)
                    if (i != k && j != k && restrict_pair[k].second == restrict_pair[i].second &&
                        restrict_pair[j].second == restrict_pair[k].first) {
                        tag = false;
                        break;
                    }
            if (tag == false) break;
        }
        if (tag == false) {
            restrict_pair.erase(restrict_pair.begin() + i);
        } else ++i;
    }


    int max_prefix_num = size * (size - 1) / 2;
    memset(restrict_last, -1, size * sizeof(int));
    memset(restrict_next, -1, max_prefix_num * sizeof(int));
    total_restrict_num = 0;
    for (const auto &p: restrict_pair) {
        // p.first must be greater than p.second
        restrict_index[total_restrict_num] = p.first;
        restrict_next[total_restrict_num] = restrict_last[p.second];
        restrict_last[p.second] = total_restrict_num;
        ++total_restrict_num;
    }
}

int Schedule::get_multiplicity() const {
    std::vector<std::vector<int> > isomorphism_vec = get_isomorphism_vec();
    return isomorphism_vec.size();
}

void Schedule::aggressive_optimize(std::vector<std::pair<int, int> > &ordered_pairs) const {
    std::vector<std::vector<int> > isomorphism_vec = get_isomorphism_vec();

    std::vector<std::vector<std::vector<int> > > permutation_groups;
    permutation_groups.clear();
    for (const std::vector<int> &v: isomorphism_vec)
        permutation_groups.push_back(calc_permutation_group(v, size));

    ordered_pairs.clear();

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (unsigned int i = 0; i < permutation_groups.size();) {
        int two_element_number = 0;
        std::pair<int, int> found_pair;
        for (const std::vector<int> &v: permutation_groups[i])
            if (v.size() == 2) {
                ++two_element_number;
                found_pair = std::pair<int, int>(v[0], v[1]);
            } else if (v.size() != 1) {
                two_element_number = -1;
                break;
            }
        if (two_element_number == 1) {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
            ordered_pairs.push_back(found_pair);
            assert(found_pair.first < found_pair.second);
        } else
            ++i;
    }

    Pattern base_dag(size);
    for (const std::pair<int, int> &pair: ordered_pairs)
        base_dag.add_ordered_edge(pair.first, pair.second);

    bool changed = true;
    while (changed && isomorphism_vec.size() != 1) {
        // use restrictions to delete other isomophism
        for (unsigned int i = 0; i < isomorphism_vec.size();) {
            Pattern test_dag(base_dag);
            const std::vector<int> &iso = isomorphism_vec[i];
            for (const std::pair<int, int> &pair: ordered_pairs)
                test_dag.add_ordered_edge(iso[pair.first], iso[pair.second]);
            if (test_dag.is_dag() == false) // is not dag means conflict
            {
                permutation_groups.erase(permutation_groups.begin() + i);
                isomorphism_vec.erase(isomorphism_vec.begin() + i);
            } else
                ++i;
        }

        changed = false;
        std::pair<int, int> found_pair;
        for (unsigned int i = 0; i < permutation_groups.size();) {
            int two_element_number = 0;
            for (const std::vector<int> &v: permutation_groups[i])
                if (v.size() == 2) {
                    ++two_element_number;
                    found_pair = std::pair<int, int>(v[0], v[1]);
                    break;
                }
            if (two_element_number >= 1) {
                permutation_groups.erase(permutation_groups.begin() + i);
                isomorphism_vec.erase(isomorphism_vec.begin() + i);
                assert(found_pair.first < found_pair.second);
                ordered_pairs.push_back(found_pair);
                base_dag.add_ordered_edge(found_pair.first, found_pair.second);
                changed = true;
                break;
            } else
                ++i;
        }
    }
    assert(isomorphism_vec.size() == 1);
}

// Schedule::aggressive_optimize(...) can only get one valid restrictions
// but in this function, we try our best to find more restrictions
// WARNING: the restrictions in ordered_pairs_vector may NOT CORRECT
void
Schedule::aggressive_optimize_get_all_pairs(std::vector<std::vector<std::pair<int, int> > > &ordered_pairs_vector) {
    std::vector<std::vector<int> > isomorphism_vec = get_isomorphism_vec();

    std::vector<std::vector<std::vector<int> > > permutation_groups;
    permutation_groups.clear();
    for (const std::vector<int> &v: isomorphism_vec)
        permutation_groups.push_back(calc_permutation_group(v, size));

    ordered_pairs_vector.clear();

    std::vector<std::pair<int, int> > ordered_pairs;
    ordered_pairs.clear();

    // delete permutation group which contains 1 permutation with 2 elements and some permutation with 1 elements,
    // and record the corresponding restriction.
    for (unsigned int i = 0; i < permutation_groups.size();) {
        int two_element_number = 0;
        std::pair<int, int> found_pair;
        for (const std::vector<int> &v: permutation_groups[i])
            if (v.size() == 2) {
                ++two_element_number;
                found_pair = std::pair<int, int>(v[0], v[1]);
            } else if (v.size() != 1) {
                two_element_number = -1;
                break;
            }
        if (two_element_number == 1) {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
            ordered_pairs.push_back(found_pair);
            assert(found_pair.first < found_pair.second);
        } else
            ++i;
    }

    Pattern base_dag(size);
    for (const std::pair<int, int> &pair: ordered_pairs)
        base_dag.add_ordered_edge(pair.first, pair.second);

    aggressive_optimize_dfs(base_dag, isomorphism_vec, permutation_groups, ordered_pairs, ordered_pairs_vector);

}

void Schedule::aggressive_optimize_dfs(Pattern base_dag, std::vector<std::vector<int> > isomorphism_vec,
                                       std::vector<std::vector<std::vector<int> > > permutation_groups,
                                       std::vector<std::pair<int, int> > ordered_pairs,
                                       std::vector<std::vector<std::pair<int, int> > > &ordered_pairs_vector) {

    for (unsigned int i = 0; i < isomorphism_vec.size();) {
        Pattern test_dag(base_dag);
        const std::vector<int> &iso = isomorphism_vec[i];
        for (const std::pair<int, int> &pair: ordered_pairs)
            test_dag.add_ordered_edge(iso[pair.first], iso[pair.second]);
        if (test_dag.is_dag() == false) // is not dag means conflict
        {
            permutation_groups.erase(permutation_groups.begin() + i);
            isomorphism_vec.erase(isomorphism_vec.begin() + i);
        } else
            ++i;
    }

    if (isomorphism_vec.size() == 1) {
        ordered_pairs_vector.push_back(ordered_pairs);
        return;
    }


    std::pair<int, int> found_pair;
    for (unsigned int i = 0; i < permutation_groups.size();) {
        int two_element_number = 0;
        for (const std::vector<int> &v: permutation_groups[i])
            if (v.size() == 2) {
                ++two_element_number;
                found_pair = std::pair<int, int>(v[0], v[1]);
                std::vector<std::vector<std::vector<int> > > next_permutation_groups = permutation_groups;
                std::vector<std::vector<int> > next_isomorphism_vec = isomorphism_vec;
                std::vector<std::pair<int, int> > next_ordered_pairs = ordered_pairs;
                Pattern next_base_dag = base_dag;

                next_permutation_groups.erase(next_permutation_groups.begin() + i);
                next_isomorphism_vec.erase(next_isomorphism_vec.begin() + i);
                assert(found_pair.first < found_pair.second);
                next_ordered_pairs.push_back(found_pair);
                next_base_dag.add_ordered_edge(found_pair.first, found_pair.second);

                aggressive_optimize_dfs(next_base_dag, next_isomorphism_vec, next_permutation_groups,
                                        next_ordered_pairs, ordered_pairs_vector);
            }
        if (two_element_number >= 1) {
            break;
        } else {
            ++i;
        }
    }

}

void Schedule::GraphZero_aggressive_optimize(std::vector<std::pair<int, int> > &ordered_pairs) const {
    std::vector<std::vector<int> > Aut;
    GraphZero_get_automorphisms(Aut);

    std::vector<std::pair<int, int> > L;
    L.clear();

    for (int v = 0; v < size; ++v) { // iterate all elements in schedule
        std::vector<std::vector<int> > stabilized_aut;
        stabilized_aut.clear();

        for (int i = 0; i < Aut.size(); ++i) {
            std::vector<int> &x = Aut[i];
            if (x[v] == v) {
                stabilized_aut.push_back(x);
            } else {
                int x1 = v, x2 = x[v];
                if (x1 > x2) {
                    int tmp = x1;
                    x1 = x2;
                    x2 = tmp;
                }
                bool tag = true;
                std::pair<int, int> cur = std::make_pair(x1, x2);
                for (int j = 0; j < L.size(); ++j)
                    if (L[j] == cur) {
                        tag = false;
                        break;
                    }
                if (tag) {
                    L.push_back(cur);
                }
            }
        }
        Aut = stabilized_aut;
    }

    ordered_pairs.clear(); // In GraphZero paper, this vector's name is 'L'

    for (int i = 0; i < L.size(); ++i) {
        bool tag = true;
        for (int j = 0; j < ordered_pairs.size(); ++j)
            if (L[i].second == ordered_pairs[j].second) {
                tag = false;
                if (L[i].first > ordered_pairs[j].first) ordered_pairs[j].first = L[i].first;
                break;
            }
        if (tag) ordered_pairs.push_back(L[i]);
    }
}

void Schedule::GraphZero_get_automorphisms(std::vector<std::vector<int> > &Aut) const {
    int p[size];
    Aut.clear();
    for (int i = 0; i < size; ++i) p[i] = i;
    do {
        bool tag = true;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j)
                if (adj_mat[INDEX(i, j, size)] != adj_mat[INDEX(p[i], p[j], size)]) {
                    tag = false;
                    break;
                }
            if (!tag) break;
        }
        if (tag) {
            std::vector<int> tmp;
            tmp.clear();
            for (int i = 0; i < size; ++i) tmp.push_back(p[i]);
            Aut.push_back(tmp);
        }
    } while (std::next_permutation(p, p + size));

}

std::vector<std::vector<int> > Schedule::get_isomorphism_vec() const {
    unsigned int pow = 1;
    for (int i = 2; i <= size; ++i)
        pow *= i;
    std::vector<std::vector<int> > vec;
    vec.clear();
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector<int> tmp_vec;
    get_full_permutation(vec, use, tmp_vec, 0);
    assert(vec.size() == pow);
    std::vector<std::vector<int> > isomorphism_vec;
    isomorphism_vec.clear();
    for (const std::vector<int> &v: vec) {
        bool flag = true;
        for (int i = 0; i < size; ++i)
            for (int j = i + 1; j < size; ++j)
                if (adj_mat[INDEX(i, j, size)] != 0)
                    if (adj_mat[INDEX(v[i], v[j], size)] == 0) // not isomorphism
                    {
                        flag = false;
                        break;
                    }
        if (flag == true)
            isomorphism_vec.push_back(v);
    }
    return isomorphism_vec;
}

void Schedule::get_full_permutation(std::vector<std::vector<int> > &vec, bool use[], std::vector<int> tmp_vec,
                                    int depth) const {
    if (depth == size) {
        vec.push_back(tmp_vec);
        return;
    }
    for (int i = 0; i < size; ++i)
        if (use[i] == false) {
            use[i] = true;
            tmp_vec.push_back(i);
            get_full_permutation(vec, use, tmp_vec, depth + 1);
            tmp_vec.pop_back();
            use[i] = false;
        }
}

std::vector<std::vector<int> > Schedule::calc_permutation_group(const std::vector<int> vec, int size) {
    bool use[size];
    for (int i = 0; i < size; ++i)
        use[i] = false;
    std::vector<std::vector<int> > res;
    res.clear();
    for (unsigned int i = 0; i < vec.size(); ++i)
        if (use[i] == false) {
            std::vector<int> tmp_vec;
            tmp_vec.clear();
            tmp_vec.push_back(i);
            use[i] = true;
            int x = vec[i];
            while (use[x] == false) {
                use[x] = true;
                tmp_vec.push_back(x);
                x = vec[x];
            }
            res.push_back(tmp_vec);
        }
    return res;
}

void Schedule::performance_modeling(int *best_order, std::vector<std::vector<int> > &candidates, int v_cnt,
                                    unsigned int e_cnt) {
    int *order;
    int *rank;

    double *p_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];

    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p;
    }

    order = new int[size];
    rank = new int[size];

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for (const std::vector<int> &vec: candidates) {
        for (int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for (int i = 1; i < size; ++i) {
            bool have_edge = false;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if (have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if (is_valid == false) continue;

        for (int i = 0; i < size; ++i) rank[order[i]] = i;
        int *cur_adj_mat;
        cur_adj_mat = new int[size * size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector<std::pair<int, int> > restricts;
        //TODO BUG!!!!!
        GraphZero_aggressive_optimize(restricts);
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());
        double *sum;
        sum = new double[restricts_size];
        for (int i = 0; i < restricts_size; ++i) sum[i] = 0;
        int *tmp;
        tmp = new int[size];
        for (int i = 0; i < size; ++i) tmp[i] = i;
        do {
            for (int i = 0; i < restricts_size; ++i)
                if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                } else break;
        } while (std::next_permutation(tmp, tmp + size));
        double total = 1;
        for (int i = 2; i <= size; ++i) total *= i;
        for (int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] / total;
        for (int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for (int i = 0; i < size; ++i) invariant_size[i].clear();
        for (int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            int cnt_backward = 0;
            for (int j = 0; j < i; ++j)
                if (cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;
            for (int j = i + 1; j < size; ++j)
                if (cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_backward;

            int c = cnt_forward;
            for (int j = i - 1; j >= 0; --j)
                if (cur_adj_mat[INDEX(j, i, size)])
                    invariant_size[j].push_back(c--);

            for (int j = 0; j < invariant_size[i].size(); ++j)
                if (invariant_size[i][j] > 1)
                    val += p_size[invariant_size[i][j] - 1] + p_size[1];
            for (int j = 0; j < restricts_size; ++j)
                if (restricts[j].second == i)
                    val *= sum[j];
            val *= p_size[cnt_forward];

        }
        if (have_best == false || val < min_val) {
            have_best = true;
            for (int i = 0; i < size; ++i)
                best_order[i] = order[i];
            min_val = val;
        }
        delete[] sum;
        delete[] tmp;
        delete[] cur_adj_mat;

    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
}

void Schedule::bug_performance_modeling(int *best_order, std::vector<std::vector<int> > &candidates, int v_cnt,
                                        unsigned int e_cnt) {
    int *order;
    int *rank;

    double *p_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];

    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p;
    }

    order = new int[size];
    rank = new int[size];

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for (const std::vector<int> &vec: candidates) {
        for (int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for (int i = 1; i < size; ++i) {
            bool have_edge = false;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if (have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if (is_valid == false) continue;

        for (int i = 0; i < size; ++i) rank[order[i]] = i;
        int *cur_adj_mat;
        cur_adj_mat = new int[size * size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector<std::vector<std::pair<int, int> > > restricts_vector;
        restricts_generate(cur_adj_mat, restricts_vector);
        for (int restricts_rank = 0; restricts_rank < restricts_vector.size(); ++restricts_rank) {
            std::vector<std::pair<int, int> > &restricts = restricts_vector[restricts_rank];
            int restricts_size = restricts.size();
            std::sort(restricts.begin(), restricts.end());
            double *sum;
            sum = new double[restricts_size];
            for (int i = 0; i < restricts_size; ++i) sum[i] = 0;
            int *tmp;
            tmp = new int[size];
            for (int i = 0; i < size; ++i) tmp[i] = i;
            do {
                for (int i = 0; i < restricts_size; ++i)
                    if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                        sum[i] += 1;
                    } else break;
            } while (std::next_permutation(tmp, tmp + size));
            double total = 1;
            for (int i = 2; i <= size; ++i) total *= i;
            for (int i = 0; i < restricts_size; ++i)
                sum[i] = sum[i] / total;
            for (int i = restricts_size - 1; i > 0; --i)
                sum[i] /= sum[i - 1];

            double val = 1;
            for (int i = 0; i < size; ++i) invariant_size[i].clear();
            for (int i = size - 1; i >= 0; --i) {
                int cnt_forward = 0;
                int cnt_backward = 0;
                for (int j = 0; j < i; ++j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_forward;
                for (int j = i + 1; j < size; ++j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_backward;

                int c = cnt_forward;
                for (int j = i - 1; j >= 0; --j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        invariant_size[j].push_back(c--);

                for (int j = 0; j < invariant_size[i].size(); ++j)
                    if (invariant_size[i][j] > 1)
                        val += p_size[invariant_size[i][j] - 1] + p_size[1];
                for (int j = 0; j < restricts_size; ++j)
                    if (restricts[j].second == i)
                        val *= sum[j];
                val *= p_size[cnt_forward];

            }
            if (have_best == false || val < min_val) {
                have_best = true;
                for (int i = 0; i < size; ++i)
                    best_order[i] = order[i];
                min_val = val;
            }
            delete[] sum;
            delete[] tmp;
        }
        delete[] cur_adj_mat;

    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
}

void Schedule::new_performance_modeling(int *best_order, std::vector<std::vector<int> > &candidates, int v_cnt,
                                        unsigned int e_cnt, long long tri_cnt) {
    int *order;
    int *rank;

    double *p_size;
    double *pp_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];
    pp_size = new double[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p0;
    }
    pp_size[0] = 1;
    for (int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i - 1] * p1;
    }

    order = new int[size];
    rank = new int[size];

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];
    for (const std::vector<int> &vec: candidates) {
        for (int i = 0; i < size; ++i)
            order[i] = vec[i];
        // check whether it is valid schedule
        bool is_valid = true;
        for (int i = 1; i < size; ++i) {
            bool have_edge = false;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if (have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if (is_valid == false) continue;

        for (int i = 0; i < size; ++i) rank[order[i]] = i;
        int *cur_adj_mat;
        cur_adj_mat = new int[size * size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector<std::vector<std::pair<int, int> > > restricts_vector;
        restricts_generate(cur_adj_mat, restricts_vector);
        for (int restricts_rank = 0; restricts_rank < restricts_vector.size(); ++restricts_rank) {
            std::vector<std::pair<int, int> > &restricts = restricts_vector[restricts_rank];
            int restricts_size = restricts.size();
            std::sort(restricts.begin(), restricts.end());
            double *sum;
            sum = new double[restricts_size];
            for (int i = 0; i < restricts_size; ++i) sum[i] = 0;
            int *tmp;
            tmp = new int[size];
            for (int i = 0; i < size; ++i) tmp[i] = i;
            do {
                for (int i = 0; i < restricts_size; ++i)
                    if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                        sum[i] += 1;
                    } else break;
            } while (std::next_permutation(tmp, tmp + size));
            double total = 1;
            for (int i = 2; i <= size; ++i) total *= i;
            for (int i = 0; i < restricts_size; ++i)
                sum[i] = sum[i] / total;
            for (int i = restricts_size - 1; i > 0; --i)
                sum[i] /= sum[i - 1];

            double val = 1;
            for (int i = 0; i < size; ++i) invariant_size[i].clear();
            for (int i = size - 1; i >= 0; --i) {
                int cnt_forward = 0;
                int cnt_backward = 0;
                for (int j = 0; j < i; ++j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_forward;
                for (int j = i + 1; j < size; ++j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        ++cnt_backward;

                int c = cnt_forward;
                for (int j = i - 1; j >= 0; --j)
                    if (cur_adj_mat[INDEX(j, i, size)])
                        invariant_size[j].push_back(c--);

                for (int j = 0; j < invariant_size[i].size(); ++j)
                    if (invariant_size[i][j] > 1)
                        val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
                val += 1;
                for (int j = 0; j < restricts_size; ++j)
                    if (restricts[j].second == i)
                        val *= sum[j];
                val *= p_size[1] * pp_size[cnt_forward - 1];

            }
            if (have_best == false || val < min_val) {
                have_best = true;
                for (int i = 0; i < size; ++i)
                    best_order[i] = order[i];
                min_val = val;
            }
            delete[] sum;
            delete[] tmp;
        }
        delete[] cur_adj_mat;

    }

    delete[] order;
    delete[] rank;
    delete[] p_size;
    delete[] pp_size;
}

void Schedule::init_in_exclusion_optimize() {
    int optimize_num = in_exclusion_optimize_num;

    assert(in_exclusion_optimize_num >= 1);

    int *id;
    id = new int[optimize_num];

    int *in_exclusion_val;
    in_exclusion_val = new int[optimize_num * 2];

    for (int n = 1; n <= optimize_num; ++n) {
        DisjointSetUnion dsu(n);
        int m = n * (n - 1) / 2;

        in_exclusion_val[2 * n - 2] = 0;
        in_exclusion_val[2 * n - 1] = 0;

        if (n == 1) {
            ++in_exclusion_val[0];
            continue;
        }

        std::pair<int, int> edge[m];
        int e_cnt = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                edge[e_cnt++] = std::make_pair(i, j);

        for (int s = 0; s < (1 << m); ++s) {
            dsu.init();
            int bit_cnt = 0;
            for (int i = 0; i < m; ++i)
                if (s & (1 << i)) {
                    ++bit_cnt;
                    dsu.merge(edge[i].first, edge[i].second);
                }
            if (dsu.get_set_size() == 1) {
                if (bit_cnt & 1) ++in_exclusion_val[2 * n - 1];
                else ++in_exclusion_val[2 * n - 2];
            }
        }
    }

    in_exclusion_optimize_group.clear();
    in_exclusion_optimize_val.clear();

    get_in_exclusion_optimize_group(0, id, 0, in_exclusion_val);

    delete[] id;
    delete[] in_exclusion_val;
}

void Schedule::get_in_exclusion_optimize_group(int depth, int *id, int id_cnt, int *in_exclusion_val) {
    if (depth == in_exclusion_optimize_num) {
        int *size = new int[id_cnt];
        for (int i = 0; i < id_cnt; ++i)
            size[i] = 0;
        for (int i = 0; i < in_exclusion_optimize_num; ++i)
            size[id[i]]++;
        int val[2];
        val[0] = in_exclusion_val[size[0] * 2 - 2];
        val[1] = in_exclusion_val[size[0] * 2 - 1];
        for (int i = 1; i < id_cnt; ++i) {
            int tmp0 = val[0];
            int tmp1 = val[1];

            val[0] = tmp0 * in_exclusion_val[size[i] * 2 - 2] + tmp1 * in_exclusion_val[size[i] * 2 - 1];
            val[1] = tmp0 * in_exclusion_val[size[i] * 2 - 1] + tmp1 * in_exclusion_val[size[i] * 2 - 2];
        }

        std::vector<std::vector<int> > group;
        group.clear();
        for (int i = 0; i < id_cnt; ++i) {
            std::vector<int> cur;
            cur.clear();
            for (int j = 0; j < in_exclusion_optimize_num; ++j)
                if (id[j] == i) cur.push_back(j);
            group.push_back(cur);
        }

        in_exclusion_optimize_group.push_back(group);
        in_exclusion_optimize_val.push_back(val[0] - val[1]);

        delete[] size;
        return;
    }

    id[depth] = id_cnt;

    get_in_exclusion_optimize_group(depth + 1, id, id_cnt + 1, in_exclusion_val);

    for (int i = 0; i < id_cnt; ++i) {
        id[depth] = i;
        get_in_exclusion_optimize_group(depth + 1, id, id_cnt, in_exclusion_val);
    }
}

void Schedule::print_schedule() const {
    printf("PATTERN_SIZE=%d\n", size);
    printf("ADJ_MATRIX=");
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%d", adj_mat[INDEX(i, j, size)]);
        }
        printf(" ");
    }
    printf("\n");
    auto tmp_pair = restrict_pair;
    std::sort(tmp_pair.begin(), tmp_pair.end(),
              [=](std::pair<int, int> &a, std::pair<int, int> &b) -> bool {
                  if (a.first == b.first)
                      return a.second < b.second;
                  else
                      return a.first < b.first;
              });

    printf("RESTRICTS=");
    for (auto restrict: tmp_pair) {
        printf("(%d,%d)", restrict.first, restrict.second);
    }
    printf("\n");
    if (in_exclusion_optimize_num > 1) {
        printf("IEP Number\t%d", in_exclusion_optimize_num);
        printf("\nIEP Redundency\t%lld", in_exclusion_optimize_redundancy);
        printf("\nIEP Val\t");
        for (auto x: in_exclusion_optimize_val) {
            printf("%d ", x);
        }
        printf("\nIEP Group\t");
        for (auto x_vec_vec: in_exclusion_optimize_group) {
            printf("[");
            int j = 0;
            for (auto x_vec: x_vec_vec){
                printf("(");
                int i = 0;
                for (auto x: x_vec) {
                    printf("%d", x);
                    if (i++ != x_vec.size() - 1) {
                        printf(" ");
                    }
                }
                printf(")");
                if (j++ != x_vec_vec.size() - 1) {
                    printf(",");
                }
            }
            printf("] | ");
        }
    }
}

void Schedule::GraphZero_performance_modeling(int *best_order, int v_cnt, unsigned int e_cnt) {
    int *order;
    int *rank;

    double *p_size;
    double *anti_p;
    p_size = new double[size];
    anti_p = new double[size];

    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    printf("fuck p %.6lf\n", p);
    p_size[0] = v_cnt;
    for (int i = 1; i < size; ++i) {
        p_size[i] = p_size[i - 1] * p;
        printf("p %d %.6lf\n", i, p_size[i]);
    }
    anti_p[0] = 1;
    for (int i = 1; i < size; ++i) {
        anti_p[i] = anti_p[i - 1] * (1 - p);
        printf("anti %d %.6lf\n", i, anti_p[i]);
    }

    order = new int[size];
    rank = new int[size];

    for (int i = 0; i < size; ++i) order[i] = i;
    double min_val;
    bool have_best = false;
    do {
        // check whether it is valid schedule
        bool is_valid = true;
        for (int i = 1; i < size; ++i) {
            bool have_edge = false;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(order[i], order[j], size)]) {
                    have_edge = true;
                    break;
                }
            if (have_edge == false) {
                is_valid = false;
                break;
            }
        }
        if (is_valid == false) continue;

        for (int i = 0; i < size; ++i) rank[order[i]] = i;
        int *cur_adj_mat;
        cur_adj_mat = new int[size * size];
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

        std::vector<std::pair<int, int> > restricts;
        GraphZero_aggressive_optimize(restricts);
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());
        double *sum;
        sum = new double[restricts_size];
        for (int i = 0; i < restricts_size; ++i) sum[i] = 0;
        int *tmp;
        tmp = new int[size];
        for (int i = 0; i < size; ++i) tmp[i] = i;
        do {
            for (int i = 0; i < restricts_size; ++i)
                if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                } else break;
        } while (std::next_permutation(tmp, tmp + size));
        double total = 1;
        for (int i = 2; i <= size; ++i) total *= i;
        for (int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] / total;
        for (int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for (int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            for (int j = 0; j < i; ++j)
                if (cur_adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;

            for (int j = 0; j < restricts_size; ++j)
                if (restricts[j].second == i)
                    val *= sum[j];
            //      val *= p_size[cnt_forward + 1] * anti_p[i - cnt_forward];
            val *= v_cnt;
            for (int j = 0; j < i - cnt_forward; ++j)
                val *= (1 - p);
            for (int j = 0; j < cnt_forward; ++j)
                val *= p;
        }
        if (have_best == false || val <= min_val) {
            have_best = true;
            for (int i = 0; i < size; ++i)
                best_order[i] = order[i];
            min_val = val;
            printf("gz upd %.10lf\n", val);
        }

        delete[] cur_adj_mat;
        delete[] sum;
        delete[] tmp;

    } while (std::next_permutation(order, order + size));

    delete[] order;
    delete[] rank;
    delete[] p_size;
    delete[] anti_p;
}

void Schedule::restrict_selection(int v_cnt, unsigned int e_cnt, long long tri_cnt,
                                  std::vector<std::vector<std::pair<int, int> > > ordered_pairs_vector,
                                  std::vector<std::pair<int, int> > &best_restricts) const {

    double *p_size;
    double *pp_size;
    int max_degree = get_max_degree();
    p_size = new double[max_degree];
    pp_size = new double[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p0;
    }

    pp_size[0] = 1;
    for (int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i - 1] * p1;
    }

    double min_val;
    bool have_best = false;
    std::vector<int> invariant_size[size];

    for (int cur_restricts = 0; cur_restricts < ordered_pairs_vector.size(); ++cur_restricts) {
        std::vector<std::pair<int, int> > &restricts = ordered_pairs_vector[cur_restricts];
        int restricts_size = restricts.size();
        std::sort(restricts.begin(), restricts.end());

        double *sum;
        sum = new double[restricts_size];
        for (int i = 0; i < restricts_size; ++i) sum[i] = 0;

        int *tmp;
        tmp = new int[size];
        for (int i = 0; i < size; ++i) tmp[i] = i;

        do {
            for (int i = 0; i < restricts_size; ++i)
                if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                    sum[i] += 1;
                } else break;
        } while (std::next_permutation(tmp, tmp + size));

        double total = 1;
        for (int i = 2; i <= size; ++i) total *= i;

        for (int i = 0; i < restricts_size; ++i)
            sum[i] = sum[i] / total;
        for (int i = restricts_size - 1; i > 0; --i)
            sum[i] /= sum[i - 1];

        double val = 1;
        for (int i = 0; i < size; ++i) invariant_size[i].clear();

        for (int i = size - 1; i >= 0; --i) {
            int cnt_forward = 0;
            int cnt_backward = 0;
            for (int j = 0; j < i; ++j)
                if (adj_mat[INDEX(j, i, size)])
                    ++cnt_forward;
            for (int j = i + 1; j < size; ++j)
                if (adj_mat[INDEX(j, i, size)])
                    ++cnt_backward;

            int c = cnt_forward;
            for (int j = i - 1; j >= 0; --j)
                if (adj_mat[INDEX(j, i, size)])
                    invariant_size[j].push_back(c--);

            for (int j = 0; j < invariant_size[i].size(); ++j)
                if (invariant_size[i][j] > 1)
                    val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
            val += 1;
            for (int j = 0; j < restricts_size; ++j)
                if (restricts[j].second == i)
                    val *= sum[j];
            val *= p_size[1] * pp_size[cnt_forward - 1];

        }
        if (have_best == false || val < min_val) {
            have_best = true;
            best_restricts = restricts;
            min_val = val;
        }

        delete[] sum;
        delete[] tmp;

    }

    delete[] p_size;
    delete[] pp_size;
    assert(have_best);
}

void Schedule::restricts_generate(const int *cur_adj_mat, std::vector<std::vector<std::pair<int, int> > > &restricts) {
    Schedule schedule(cur_adj_mat, get_size());
    schedule.aggressive_optimize_get_all_pairs(restricts);
    int size = schedule.get_size();
    Graph *complete;
    DataLoader *D = new DataLoader();
    assert(D->load_complete(complete, size + 1));
    int thread_num = 1;
    long long ans = complete->pattern_matching(schedule, thread_num) / schedule.get_multiplicity();
    for (int i = 0; i < restricts.size();) {
        Schedule cur_schedule(schedule.get_adj_mat_ptr(), schedule.get_size());
        cur_schedule.add_restrict(restricts[i]);
        long long cur_ans = complete->pattern_matching(cur_schedule, thread_num);
        if (cur_ans != ans) {
            restricts.erase(restricts.begin() + i);
        } else {
            ++i;
        }
    }

    delete complete;
    delete D;
}

int Schedule::get_vec_optimize_num(const std::vector<int> &vec) {
    bool is_valid = true;
    for (int i = 1; i < size; ++i) {
        bool have_edge = false;
        for (int j = 0; j < i; ++j)
            if (adj_mat[INDEX(vec[i], vec[j], size)]) {
                have_edge = true;
                break;
            }
        if (have_edge == false) {
            is_valid = false;
            break;
        }
    }
    if (!is_valid) return -1;

    for (int k = 2; k <= size; ++k) {
        bool flag = true;
        for (int i = size - k + 1; i < size; ++i)
            if (adj_mat[INDEX(vec[size - k], vec[i], size)]) {
                flag = false;
                break;
            }
        if (flag == false) return k - 1;
    }
    assert(0);
    return -1;
}

double
Schedule::our_estimate_schedule_restrict(const std::vector<int> &order, const std::vector<std::pair<int, int> > &pairs,
                                         uint64_t v_cnt, uint64_t e_cnt, uint64_t tri_cnt) {
    int max_degree = get_max_degree();

    double p_size[max_degree];
    double pp_size[max_degree];

    double p0 = e_cnt * 1.0 / v_cnt / v_cnt;
    double p1 = tri_cnt * 1.0 * v_cnt / e_cnt / e_cnt;

    p_size[0] = v_cnt;
//    printf("0th loop p %f\n", p_size[0]);

    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p0;
//        printf("%dth loop p %f\n", i, p_size[i]);
    }
    pp_size[0] = 1;
    for (int i = 1; i < max_degree; ++i) {
        pp_size[i] = pp_size[i - 1] * p1;
//        printf("%dth loop pp %f\n", i, pp_size[i]);
    }

    int rank[size];
    for (int i = 0; i < size; ++i) rank[order[i]] = i;

    int *cur_adj_mat;
    cur_adj_mat = new int[size * size];
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int, int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i) sum[i] = 0;

    int tmp[size];
    for (int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i)
            if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            } else break;
    } while (std::next_permutation(tmp, tmp + size));

    double total = 1;
    for (int i = 2; i <= size; ++i) total *= i;
    for (int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] / total;
    for (int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for (int i = 0; i < size; ++i) invariant_size[i].clear();

    double val = 1;
    for (int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for (int j = 0; j < i; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for (int j = i + 1; j < size; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for (int j = i - 1; j >= 0; --j)
            if (cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for (int j = 0; j < invariant_size[i].size(); ++j)
            if (invariant_size[i][j] > 1)
                val += p_size[1] * pp_size[invariant_size[i][j] - 2] + p_size[1];
        val += 1;
        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *= sum[j];
        if (i) {
            val *= p_size[1] * pp_size[cnt_forward - 1];
        } else {
            val *= p_size[0];
        }
    }
    delete[] cur_adj_mat;
    return val;
}

double Schedule::GraphZero_estimate_schedule_restrict(const std::vector<int> &order,
                                                      const std::vector<std::pair<int, int> > &pairs, uint64_t v_cnt,
                                                      uint64_t e_cnt) {
    int max_degree = get_max_degree();

    double p_size[max_degree];
    double p = e_cnt * 1.0 / v_cnt / v_cnt;

    p_size[0] = v_cnt;
    for (int i = 1; i < max_degree; ++i) {
        p_size[i] = p_size[i - 1] * p;
    }

    int rank[size];
    for (int i = 0; i < size; ++i) rank[order[i]] = i;

    int *cur_adj_mat;
    cur_adj_mat = new int[size * size];
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int, int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i) sum[i] = 0;

    int tmp[size];
    for (int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i) {
            if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            } else break;
        }
    } while (std::next_permutation(tmp, tmp + size));

    double total = 1;
    for (int i = 2; i <= size; ++i) total *= i;

    for (int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] / total;
    for (int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    std::vector<int> invariant_size[size];
    for (int i = 0; i < size; ++i) invariant_size[i].clear();

    double val = 1;
    for (int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        int cnt_backward = 0;
        for (int j = 0; j < i; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;
        for (int j = i + 1; j < size; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_backward;

        int c = cnt_forward;
        for (int j = i - 1; j >= 0; --j)
            if (cur_adj_mat[INDEX(j, i, size)])
                invariant_size[j].push_back(c--);

        for (int j = 0; j < invariant_size[i].size(); ++j)
            if (invariant_size[i][j] > 1)
                val += p_size[invariant_size[i][j] - 1] + p_size[1];
        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *= sum[j];
        val *= p_size[cnt_forward];

    }

    delete[] cur_adj_mat;

    return val;
}

double Schedule::Naive_estimate_schedule_restrict(const std::vector<int> &order,
                                                  const std::vector<std::pair<int, int> > &pairs, int v_cnt,
                                                  unsigned int e_cnt) {

    double p = e_cnt * 2.0 / v_cnt / v_cnt;

    int rank[size];
    for (int i = 0; i < size; ++i) rank[order[i]] = i;

    int *cur_adj_mat;
    cur_adj_mat = new int[size * size];
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            cur_adj_mat[INDEX(rank[i], rank[j], size)] = adj_mat[INDEX(i, j, size)];

    std::vector<std::pair<int, int> > restricts = pairs;
    int restricts_size = restricts.size();
    std::sort(restricts.begin(), restricts.end());

    double sum[restricts_size];
    for (int i = 0; i < restricts_size; ++i) sum[i] = 0;
    int tmp[size];

    for (int i = 0; i < size; ++i) tmp[i] = i;
    do {
        for (int i = 0; i < restricts_size; ++i)
            if (tmp[restricts[i].first] > tmp[restricts[i].second]) {
                sum[i] += 1;
            } else break;
    } while (std::next_permutation(tmp, tmp + size));

    double total = 1;
    for (int i = 2; i <= size; ++i) total *= i;

    for (int i = 0; i < restricts_size; ++i)
        sum[i] = sum[i] / total;
    for (int i = restricts_size - 1; i > 0; --i)
        sum[i] /= sum[i - 1];

    double val = 1;
    for (int i = size - 1; i >= 0; --i) {
        int cnt_forward = 0;
        for (int j = 0; j < i; ++j)
            if (cur_adj_mat[INDEX(j, i, size)])
                ++cnt_forward;

        for (int j = 0; j < restricts_size; ++j)
            if (restricts[j].second == i)
                val *= sum[j];
        val *= v_cnt;
        for (int j = 0; j < i - cnt_forward; ++j)
            val *= (1 - p);
        for (int j = 0; j < cnt_forward; ++j)
            val *= p;
    }

    delete[] cur_adj_mat;

    return val;
}

void Schedule::remove_invalid_permutation(std::vector<std::vector<int> > &candidate_permutations) {
    for (unsigned int i = 0; i < candidate_permutations.size();) {
        const auto &vec = candidate_permutations[i];
        bool tag = true;
        for (int x = 1; x < size; ++x) {
            bool have_edge = false;
            for (int y = 0; y < x; ++y)
                if (adj_mat[INDEX(vec[x], vec[y], size)]) {
                    have_edge = true;
                    break;
                }
            if (!have_edge) {
                tag = false;
                break;
            }
        }
        if (tag) {
            ++i;
        } else {
            candidate_permutations.erase(candidate_permutations.begin() + i);
        }
    }
}

int Schedule::get_in_exclusion_optimize_num_when_not_optimize() {
    std::vector<int> I;
    for (int i = 0; i < size; ++i) I.push_back(i);
    return get_vec_optimize_num(I);
}

void Schedule::set_in_exclusion_optimize_redundancy() {
    int tmp = get_in_exclusion_optimize_num();
    if (tmp <= 1) {
        in_exclusion_optimize_redundancy = 1;
    } else {
        Graph *complete;
        DataLoader *D = new DataLoader();
        assert(D->load_complete(complete, get_size()));
        delete D;
        in_exclusion_optimize_redundancy = 1;
        long long ans = complete->pattern_matching(*this, 1);
        set_in_exclusion_optimize_num(0);
        long long true_ans = complete->pattern_matching(*this, 1);
        set_in_exclusion_optimize_num(tmp);
        delete complete;
        in_exclusion_optimize_redundancy = ans / true_ans;
    }
}

