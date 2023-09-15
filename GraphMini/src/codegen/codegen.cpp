//
// Created by ubuntu on 1/1/23.
//


#include "ir.h"
#include "codegen.h"
#include "logging.h"
#include "timer.h"

#include <cmath>
#include <fmt/format.h>
#include <sstream>
#include <algorithm>
#include <utility>
#include "../../dependency/GraphPi/include/schedule.h"

namespace minigraph {
    bool EnableProfling = false;
    CodeGenConfig CurConfig;

    inline int VEC_INDEX(int i, int j, int p_size) { return i * p_size + j; };

    inline std::string gen_indent(int dep) {
        return std::string(dep + 4, '\t');
    };

    inline std::string gen_indent_tbb(int dep) {
        return std::string(dep + 4, '\t');
    };

    std::string restricts_to_str(std::vector<std::pair<int, int>> &pair_vec, int p_size) {
        std::string mat;
        mat.resize(p_size * p_size, '0');

        for (auto pair: pair_vec) {
            mat.at(VEC_INDEX(pair.second, pair.first, p_size)) = '1';
        }
        for (int v1 = 0; v1 < p_size; v1++) {
            std::vector<int> greater_than_v1;
            for (int v1_plus = 0; v1_plus < p_size; ++v1_plus) {
                if (mat.at(VEC_INDEX(v1, v1_plus, p_size)) == '1') {
                    greater_than_v1.push_back(v1_plus);
                }
            }

            for (auto &pair: pair_vec) {
                if (pair.first == v1) {
                    int less_than_v1 = pair.second;
                    for (auto v1_plus: greater_than_v1) {
                        mat.at(VEC_INDEX(less_than_v1, v1_plus, p_size)) = '1';
                    }
                }
            }
        }

        return mat;
    }

    int get_pattern_size(const std::string &adj_mat) {
        int pattern_size = (int) sqrt(adj_mat.size());
        if (pattern_size * pattern_size != adj_mat.size() || pattern_size < 3)
            LOG(FATAL) << "Invalid adj matrix size, it must be a square number";
        return pattern_size;
    }

    EdgeIR ToEdgeIR(const std::string &adj_mat, int vid) {
        if (vid == 0) return {};
        int p_size = get_pattern_size(adj_mat);
        EdgeIR out;
        for (int j = 0; j < p_size; j++) {
            out[j] = adj_mat.at(VEC_INDEX(vid, j, p_size)) == '1';
        }
        return out;
    }

    EdgeRestrictIR ToRestrictIR(const std::string &res_mat, int vid) {
        return ToEdgeIR(res_mat, vid);
    };

    PlanIR create_plan(const std::string &_adj_mat, CodeGenConfig config, MetaData meta) {
        Timer t;
        Schedule sc{};
        int p_size = get_pattern_size(_adj_mat);
        sc.get_schedule(_adj_mat.c_str(), p_size, meta.num_vertex, meta.num_edge, meta.num_triangle);
        std::string adj_mat = sc.get_adj_mat_str();
        std::string res_mat = restricts_to_str(sc.restrict_pair, p_size);
        LOG(MSG) << "SCHEDULING_TIME(s)=" << t.Passed();
        t.Reset();
        PlanIR out;
        VertexSetIR::adjMatType = config.adjMatType;
        int max_dep = p_size - 1;
        std::vector<EdgeIR> edge_ir_vec(p_size);
        std::vector<EdgeRestrictIR> res_ir_vec(p_size);

        std::vector<std::vector<VertexSetIR>> set_ops(max_dep);
        std::vector<VertexSetIR> iter_set(max_dep);

        for (int vid = 1; vid < p_size; vid++) {
            edge_ir_vec.at(vid) = ToEdgeIR(adj_mat, vid);
            res_ir_vec.at(vid) = ToRestrictIR(res_mat, vid);
        }

        for (int dep = 0; dep < max_dep; dep++) {
            for (int vid = dep + 1; vid < p_size; vid++) {
                const EdgeIR &e = edge_ir_vec.at(vid);
                const EdgeRestrictIR &r = res_ir_vec.at(vid);
                VertexSetIR v_ir = VertexSetIR(e, r, dep);
                if (v_ir.edge_num() > 0)
                    set_ops.at(dep).push_back(v_ir);
            }

            // next vertex set to iterate through
            int iter_vid = dep + 1;
            const EdgeIR &e = edge_ir_vec.at(iter_vid);
            const EdgeRestrictIR &r = res_ir_vec.at(iter_vid);
            VertexSetIR v_iter_ir = VertexSetIR(e, r, dep);
            CHECK(v_iter_ir.edge_num() > 0)
                << "Invalid schedule: The set of vertex to iterate potentially contains all vertices in the graph";
            iter_set.at(dep) = v_iter_ir;
        }
        // remove duplicated vertexIR
        int next_id = 0;
        for (auto &ops: set_ops) {
            std::sort(ops.begin(), ops.end());
            ops.erase(std::unique(ops.begin(), ops.end()), ops.end());
            for (auto &op: ops) {
                op.id = next_id++;
            }
        }

        for (int dep = 0; dep < max_dep; dep++) {
            const auto &ops = set_ops.at(dep);
            VertexSetIR iter_vs = *std::find(ops.begin(), ops.end(), iter_set.at(dep));
            iter_set.at(dep) = iter_vs;
            CHECK(iter_vs.has_id()) << "Invalid VertexSetIR (id = -1)";
        }
        out.p_size = p_size;
        out.set_ops = set_ops;
        out.iter_set = iter_set;
        // iep optimization

        if (config.adjMatType == AdjMatType::EdgeInducedIEP) {
            out.iep_num = sc.get_in_exclusion_optimize_num();
            out.iep_depth = p_size - out.iep_num - 1;
            out.iep_groups = sc.in_exclusion_optimize_group;
            out.iep_vals = sc.in_exclusion_optimize_val;
            out.iep_redundancy = sc.get_in_exclusion_optimize_redundancy();

            if (out.iep_num > 1) {
                for (int vset_id = 0; vset_id < out.iep_num; ++vset_id) {
                    int iter_depth = vset_id + out.iep_depth;
                    std::optional<VertexSetIR> vset;
                    VertexSetIR child = out.iter_set.at(iter_depth);
                    for (VertexSetIR parent: out.set_ops.at(out.iep_depth)) {
                        if (parent.same_iep_computation(child))
                            vset = parent;
                    }
                    assert(vset.has_value());
                    out.iep_set.push_back(vset.value());
                }
            }
        }
        // LOG(MSG) << "MiniGraph" << t.Passed() << "s";
        out.meta = meta;
        return out;
    };

    PlanIR create_plan_mg(const PlanIR &_plan, const CodeGenConfig &_config) {
        PlanIR out = _plan;
        int max_dep = std::min(out.p_size - 2, out.p_size - out.iep_num - 1);
        out.mg_ops.resize(out.p_size);
        for (int comp_dep = max_dep; comp_dep >= 2; --comp_dep) {
            const std::vector<VertexSetIR> &vset_vec = out.set_ops.at(comp_dep);
            const VertexSetIR &iter = out.iter_set.at(comp_dep - 1);
            for (const VertexSetIR &vset: vset_vec) {
                if (VertexSetIR::adjMatType != AdjMatType::VertexInduced && vset.is_edge(comp_dep) == false) {
                    continue; // EdgeInduced matching does not require set subtraction
                }

                for (int prune_dep = 0; prune_dep <= comp_dep - 2; prune_dep++) {
                    std::optional<VertexSetIR> intersect = out.get_parent_vset(vset, prune_dep);
                    std::optional<VertexSetIR> vertices = out.get_parent_vset(iter, prune_dep);
                    if (intersect.has_value() && vertices.has_value()) {
                        MiniGraphIR mg_op(vertices.value(), intersect.value());
                        if (_config.pruningType == PruningType::Static) {
                            VertexSetIR next_vertices = out.iter_set.at(prune_dep); // next iteration
                            for (const auto &next_intersect: out.set_ops.at(prune_dep + 1)) {
                                if (mg_op.computed(next_vertices, next_intersect)) {
                                    out.mg_ops.at(prune_dep).push_back(mg_op);
                                    break;
                                }
                            }
                        } else {
                            out.mg_ops.at(prune_dep).push_back(mg_op);
                        }
                    } // find a valid minigraph
                }
            }
        }
        // remove duplicates
        int next_id = 0;
        for (auto &mg_op: out.mg_ops) {
            std::sort(mg_op.begin(), mg_op.end());
            mg_op.erase(std::unique(mg_op.begin(), mg_op.end()), mg_op.end());
            for (auto &mg: mg_op) {
                mg.id = next_id++;
            }
        }
        std::vector<bool> IndeedUsed(next_id, false);
        out.mg_used.resize(out.p_size);
        out.mg_bounded.clear();
        out.mg_bounded.resize(next_id, true);
        for (int dep = 1; dep <= max_dep; dep++) {
            const auto &ops = out.set_ops.at(dep);
            std::vector<MiniGraphIR> mg_used_vec;
            for (const VertexSetIR &op: ops) {
                std::optional<MiniGraphIR> mg_used = out.get_parent_mg(op);
                if (mg_used.has_value()) {
                    IndeedUsed.at(mg_used->id) = true;
                    bool to_add = true;
                    for (const auto &mg: mg_used_vec) {
                        if (mg.id == mg_used->id) to_add = false;
                    }
                    if (to_add) {
                        mg_used_vec.push_back(mg_used.value());
                        out.mg_bounded.at(mg_used->id) =
                                out.mg_bounded.at(mg_used->id) && op.is_restricted(op.loop_depth());
                    }
                }
            }
            out.mg_used.at(dep) = mg_used_vec;
        }

        for (int dep = max_dep; dep >= 0; dep--) {
            for (const MiniGraphIR &child: out.mg_ops.at(dep)) {
                if (IndeedUsed.at(child.id)) {
                    std::optional<MiniGraphIR> parent = out.get_parent_mg(child);
                    if (parent.has_value()) {
                        IndeedUsed.at(parent->id) = true;
                        out.mg_bounded.at(parent->id) = out.mg_bounded.at(parent->id) && out.mg_bounded.at(child.id);
                    }
                }
            }
        }
        for (auto &mg_op: out.mg_ops) {
            for (auto itr = mg_op.begin(); itr != mg_op.end(); itr++) {
                if (!IndeedUsed.at(itr->id)) {
                    mg_op.erase(itr);
                }
            }
        }
        return out;
    };


    std::string gen_code_read_adj(const PlanIR &plan, int dep) {
        std::string out;
        if (EnableProfling) {
            out += fmt::format("ctx.profiler->set_cur_loop({dep});\n", fmt::arg("dep", dep));
            out += gen_indent(dep);
        }
        bool NoAdjNeeded = true;
        if (CurConfig.pruningType == PruningType::None) {
            NoAdjNeeded = false;
        } else {
            for (const auto &op: plan.set_ops.at(dep)) {
                auto mg = plan.get_parent_mg(op);
                NoAdjNeeded = NoAdjNeeded && mg.has_value();
            }
        }
        if (dep > 0) {
            const VertexSetIR &iter = plan.iter_set.at(dep - 1);
            out += fmt::format("const IdType i{dep}_id = s{iter_id}[i{dep}_idx];\n",
                                                 fmt::arg("dep", dep),
                                                 fmt::arg("iter_id", iter.id));
        }

        if (!NoAdjNeeded) {
            if (dep > 0) out += gen_indent(dep);
            out += fmt::format("VertexSet i{dep}_adj = graph->N(i{dep}_id);\n",
                                                 fmt::arg("dep", dep));
        }
        return out;
    }

    std::string gen_code_iter(const PlanIR &plan, int dep) {
        if (dep >= plan.p_size - 2) {
            return "";
        } else {
            const auto &iter_set = plan.iter_set.at(dep);
            return fmt::format(
                    "for (size_t i{dep}_idx = 0; i{dep}_idx < s{iter_id}.size(); i{dep}_idx++) {left} // loop-{dep} begin\n",
                    fmt::arg("left", "{"),
                    fmt::arg("iter_id", iter_set.id),
                    fmt::arg("dep", dep + 1));
        }
    }

    std::string gen_code_op(const PlanIR &plan, const VertexSetIR &op) {
        std::string out;
        auto parent = plan.get_parent_vset(op);
        int dep = op.loop_depth();
        std::string upper_bound;
        if (op.is_restricted(op.loop_depth())) {
            upper_bound = fmt::format(", i{}_adj.vid()", dep);
        }
        if (parent.has_value()) {
            CHECK(parent->loop_depth() + 1 >= op.loop_depth()) << "Not optimal parent";
            // generate code from prefix
            if (parent->loop_depth() == op.loop_depth()) {
                CHECK(op.is_restricted(op.loop_depth())) << "\nLogic error (op is not restricted at loop_depth)\nOP:\n"
                                                         << op << "\nParent:\n" << parent.value();
                out += fmt::format("VertexSet s{op_id} = s{parent_id}.bounded(i{dep}_id);\n",
                                   fmt::arg("op_id", op.id),
                                   fmt::arg("parent_id", parent->id),
                                   fmt::arg("dep", dep));

            } else if (parent->loop_depth() == op.loop_depth() - 1) {

                if (op.is_edge(op.loop_depth())) {
                    // Both VertexInduced and EdgeInduced: intersection
                    if (!plan.is_last_op(op)) {
                        // intersect and return vertex set
                        out += fmt::format(
                                "VertexSet s{op_id} = s{parent_id}.intersect(i{dep}_adj{upper_bound});\n",
                                fmt::arg("op_id", op.id),
                                fmt::arg("parent_id", parent->id),
                                fmt::arg("dep", dep),
                                fmt::arg("upper_bound", upper_bound));
                    } else {
                        // intersect and return counter
                        out += fmt::format("counter += s{parent_id}.intersect_cnt(i{dep}_adj{upper_bound});\n",
                                           fmt::arg("parent_id", parent->id),
                                           fmt::arg("dep", dep),
                                           fmt::arg("upper_bound", upper_bound));
                    }
                } else {
                    if (VertexSetIR::adjMatType == minigraph::AdjMatType::VertexInduced) {
                        // VertexInduced: subtraction
                        if (!plan.is_last_op(op)) {
                            // last op: intersect and return vertex set
                            out += fmt::format(
                                    "VertexSet s{op_id} = s{parent_id}.subtract(i{dep}_adj{upper_bound});\n",
                                    fmt::arg("op_id", op.id),
                                    fmt::arg("parent_id", parent->id),
                                    fmt::arg("dep", dep),
                                    fmt::arg("upper_bound", upper_bound));
                        } else {
                            // intersect and return counter
                            out += fmt::format("counter += s{parent_id}.subtract_cnt(i{dep}_adj{upper_bound});\n",
                                               fmt::arg("parent_id", parent->id),
                                               fmt::arg("dep", dep),
                                               fmt::arg("upper_bound", upper_bound));
                        }
                    } else {
                        // EdgeInduced: remove v_iter
                        if (!plan.is_last_op(op)) {
                            if (op.is_restricted(op.loop_depth())) {
                                out += fmt::format("VertexSet s{op_id} = s{parent_id}.bounded(i{dep}_adj.vid());\n",
                                                   fmt::arg("op_id", op.id),
                                                   fmt::arg("parent_id", parent->id),
                                                   fmt::arg("dep", dep));
                            } else {
                                out += fmt::format("VertexSet s{op_id} = s{parent_id}.remove(i{dep}_adj.vid());\n",
                                                   fmt::arg("op_id", op.id),
                                                   fmt::arg("parent_id", parent->id),
                                                   fmt::arg("dep", dep));
                            }
                        } else {
                            if (op.is_restricted(op.loop_depth())) {
                                out += fmt::format("counter += s{parent_id}.bounded_cnt(i{dep}_adj.vid());\n",
                                                   fmt::arg("op_id", op.id),
                                                   fmt::arg("parent_id", parent->id),
                                                   fmt::arg("dep", dep));
                            } else {
                                out += fmt::format("counter += s{parent_id}.remove_cnt(i{dep}_adj.vid());\n",
                                                   fmt::arg("op_id", op.id),
                                                   fmt::arg("parent_id", parent->id),
                                                   fmt::arg("dep", dep));
                            }
                        }
                    }

                }
            }
        } else {
            // no parent
            CHECK(op.edge_num() == 1 && op.is_edge(op.loop_depth()))
                << "\nLogic error: VertexSetIR should have one parent but get none\n" << op;

            if (op.is_restricted(op.loop_depth())) {
                out += fmt::format("VertexSet s{op_id} = i{dep}_adj.bounded(i{dep}_id)",
                                   fmt::arg("op_id", op.id),
                                   fmt::arg("dep", dep));
            } else {
                out += fmt::format("VertexSet s{op_id} = i{dep}_adj",
                                   fmt::arg("op_id", op.id),
                                   fmt::arg("dep", dep));
            }

            for (int subtract_id = 0; subtract_id < dep; subtract_id++) {
                if (VertexSetIR::adjMatType == minigraph::AdjMatType::VertexInduced) {
                    std::string subtract_bound;
                    if (op.is_restricted(subtract_id)) {
                        subtract_bound = fmt::format(", i{}_adj.vid()", subtract_id);
                    }
//                    out += fmt::format("s{op_id} = s{op_id}.subtract(i{subtract_id}_adj{upper_bound});\n",
//                                   fmt::arg("op_id", op.id),
//                                   fmt::arg("iter_id", dep),
//                                   fmt::arg("subtract_id", subtract_id),
//                                   fmt::arg("upper_bound", subtract_bound));
                    out += fmt::format(".subtract(i{subtract_id}_adj{upper_bound})",
                                       fmt::arg("op_id", op.id),
                                       fmt::arg("iter_id", dep),
                                       fmt::arg("subtract_id", subtract_id),
                                       fmt::arg("upper_bound", subtract_bound));

                } else { // EdgeInduced
                    if (op.is_restricted(subtract_id)) {
//                        out += fmt::format("s{op_id} = s{op_id}.bounded(i{subtract_id}_adj.vid());\n",
//                                           fmt::arg("op_id", op.id),
//                                           fmt::arg("iter_id", dep),
//                                           fmt::arg("subtract_id", subtract_id));
                        out += fmt::format(".bounded(i{subtract_id}_adj.vid())",
                                           fmt::arg("op_id", op.id),
                                           fmt::arg("iter_id", dep),
                                           fmt::arg("subtract_id", subtract_id));
                    } else {
//                        out += fmt::format("s{op_id} = s{op_id}.remove(i{subtract_id}_adj.vid());\n",
//                                           fmt::arg("op_id", op.id),
//                                           fmt::arg("iter_id", dep),
//                                           fmt::arg("subtract_id", subtract_id));
                        out += fmt::format(".remove(i{subtract_id}_adj.vid())",
                                           fmt::arg("op_id", op.id),
                                           fmt::arg("iter_id", dep),
                                           fmt::arg("subtract_id", subtract_id));
                    }
                }
            }
            out += ";\n";
            out += gen_indent(dep) + fmt::format("if (s{op_id}.size() == 0) continue;\n", fmt::arg("op_id", op.id));
            if (plan.is_last_op(op)) {
                out += fmt::format("counter += s{op_id}.size();\n", fmt::arg("op_id", op.id));
            }
        }

        return out;
    }

    bool mg_should_eager(const PlanIR& plan, const MiniGraphIR& mg){
        // for (int loop_dep = mg.loop_depth() + 1; loop_dep < plan.mg_ops.size(); loop_dep++){
        //     for (auto& cmg: plan.mg_ops.at(loop_dep)) {
        //         if (mg.is_superset_of(cmg)) return true;
        //     }
        // }

        // if (plan.get_parent_mg(mg).has_value()) return true;

        VertexSetIR next_vertices = plan.iter_set.at(mg.loop_depth()); // next iteration
        if (!(next_vertices == mg.m_vertices)) return false;
        for (const auto &next_intersect: plan.set_ops.at(mg.loop_depth() + 1)) {
            if (mg.computed(next_vertices, next_intersect)) {
                return true;
            }
        }
        return false;
    }

    std::string gen_mg_type(const PlanIR& plan, const MiniGraphIR &mg) {
        return mg_should_eager(plan, mg) ? "MiniGraphEager" : "MiniGraphType";
    }

    std::string gen_code_mg_init(const PlanIR &plan, const MiniGraphIR &mg) {
//        const VertexSetIR &iter = plan.iter_set.at(mg.loop_depth());
        std::string mgType = gen_mg_type(plan, mg);
        return fmt::format("{mg_type} m{mg_id}({_bounded},{_par});\n",
                           fmt::arg("mg_id", mg.id),
                           fmt::arg("mg_type", mgType),
                           fmt::arg("_bounded", plan.is_bounded(mg)),
                           fmt::arg("_par", plan.is_par(mg)));
    }

    std::string gen_code_mg_adj(const PlanIR &plan, int dep, int indent_dep = -1) {
        if (dep == 0) return "";
        std::string indent = (indent_dep == -1) ? gen_indent(dep) : gen_indent_tbb(dep);
        const VertexSetIR &iter = plan.iter_set.at(dep - 1);
        std::string out;
        for (const auto &mg: plan.mg_used.at(dep)) {
            auto &vertices = mg.m_vertices;
            bool same_address = true;
            if (iter.loop_depth() == vertices.loop_depth()) {
                for (int i = 0; i <= iter.loop_depth(); i++) {
                    if (iter.is_edge(i) != vertices.is_edge(i)) same_address = false;
                }
            } else {
                same_address = false;
            }

            if (same_address) {
                out += indent;
                out += fmt::format("VertexSet m{mg_id}_adj = m{mg_id}.N(i{dep}_idx);\n",
                                   fmt::arg("mg_id", mg.id),
                                   fmt::arg("dep", dep));
            } else {
                std::string v_idx = fmt::format("m{mg_id}_s{iter_id}[i{dep}_idx]",
                                                fmt::arg("mg_id", mg.id),
                                                fmt::arg("iter_id", iter.id),
                                                fmt::arg("dep", dep));
                out += indent;
                out += fmt::format("VertexSet m{mg_id}_adj = m{mg_id}.N({v_idx});\n",
                                   fmt::arg("mg_id", mg.id),
                                   fmt::arg("v_idx", v_idx));
            }
        }
        return out;
    }

    bool skip_build_indices(const PlanIR &plan, const MiniGraphIR &mg, const VertexSetIR &iter) {
        auto &vertices = mg.m_vertices;
        if (iter.loop_depth() == vertices.loop_depth()) {
            bool same_address = true;
            for (int i = 0; i <= iter.loop_depth(); i++) {
                if (iter.is_edge(i) != vertices.is_edge(i)) same_address = false;
            }
            if (same_address) return true;
        }
        return false;
    };

    std::string gen_code_mg_indice(const PlanIR &plan, const MiniGraphIR &mg, int dep) {
        const VertexSetIR &iter = plan.iter_set.at(dep);
        if (skip_build_indices(plan, mg, iter))
            return fmt::format("//skip building indices for m{mg_id} because they can be obtained directly\n",
                               fmt::arg("mg_id", mg.id)); // skip
        else
            return fmt::format("auto m{mg_id}_s{iter_id} = m{mg_id}.indices(s{iter_id});\n",
                               fmt::arg("mg_id", mg.id),
                               fmt::arg("iter_id", iter.id));
    };

    std::string gen_code_mg_est_visits(const PlanIR &plan, const MiniGraphIR &mg, int iter_dep) {
        int iter_id = plan.iter_set.at(mg.loop_depth()).id;
        std::string out = fmt::format("s{}.size()", iter_id);
        for (int dep = mg.loop_depth() + 1; dep < iter_dep; dep++) {
            const VertexSetIR &cur_iter = plan.iter_set.at(dep);
            auto parent = plan.get_parent_vset(cur_iter, mg.loop_depth());
            if (parent.has_value()) {
                double p1 = 1.0 * plan.meta.num_edge / plan.meta.num_vertex / plan.meta.num_vertex;
                double p2 = 1.0 * plan.meta.num_triangle * 6 * plan.meta.num_vertex / plan.meta.num_edge / plan.meta.num_edge;
                double rate = 1.0;
                for (int adj_dep = mg.loop_depth(); adj_dep < cur_iter.loop_depth(); adj_dep++){
                    auto &adj_iter = plan.iter_set.at(adj_dep);
                    if (adj_iter.share_at_least_one_parent_node(parent.value())) {
                        rate *= p2;
                    } else {
                        rate *= p1;
                    }
                }

                out += fmt::format(" * s{par_id}.size() * {rate}",
                                   fmt::arg("par_id", parent->id),
                                   fmt::arg("rate", rate));
            } else {
                double avg_deg = 1.0 * plan.meta.num_edge / plan.meta.num_vertex;
                out += fmt::format(" * {}", avg_deg);
            }
        }
        return out;
    }

    std::string gen_code_mg_build(const PlanIR &plan, const MiniGraphIR &mg) {
        const VertexSetIR &iter = plan.iter_set.at(mg.loop_depth());
        int iter_id = iter.id;
        std::optional<MiniGraphIR> parent_mg = plan.get_parent_mg(mg);
        std::string out;
        bool should_build_mg_eagerly = mg_should_eager(plan, mg);

        if (CurConfig.pruningType == PruningType::CostModel && !should_build_mg_eagerly) {
                std::string factor = fmt::format("double m{mg_id}_factor = 0;\n", fmt::arg("mg_id", mg.id));
                int max_dep = std::min(plan.p_size - 2, plan.p_size - plan.iep_num - 1);

                int total_reuse = 0;
                for (int dep = mg.loop_depth() + 2; dep <= max_dep; dep++) {
//                    const VertexSetIR &cur_iter = plan.iter_set.at(dep-1);
                    int multiplier = 0;
                    for (const auto &op: plan.set_ops.at(dep)) {
                        auto parent_mg = plan.get_parent_mg(op);
                        if (parent_mg.has_value()) {
                            if (mg.is_superset_of(parent_mg.value()) || mg == parent_mg.value()) multiplier++;
                        }
                    }
                    if (multiplier > 0) {
                        std::string est_visits = gen_code_mg_est_visits(plan, mg, dep);
                        factor += gen_indent(mg.loop_depth()) + fmt::format("m{mg_id}_factor += {est_visits} * {multiplier};\n",
                                                                            fmt::arg("mg_id", mg.id), fmt::arg("est_visits", est_visits),
                                                                            fmt::arg("multiplier", multiplier));
                    }
                    total_reuse += multiplier;
                };
                // assert(total_reuse > 0);
                out += factor;
                out += gen_indent(mg.loop_depth());
                out += fmt::format("m{mg_id}.set_reuse_multiplier(m{mg_id}_factor); ",
                                   fmt::arg("mg_id", mg.id));
        }

        if (parent_mg.has_value()) {
            out += fmt::format("m{mg_id}.build(&m{parent_id}, s{vset_id}, s{vint_id}, s{iter_id});\n",
                               fmt::arg("parent_id", parent_mg->id),
                               fmt::arg("mg_id", mg.id), fmt::arg("vset_id", mg.vset_id()),
                               fmt::arg("vint_id", mg.vint_id()), fmt::arg("iter_id", iter_id));
        } else {
            out += fmt::format("m{mg_id}.build(s{vset_id}, s{vint_id}, s{iter_id});\n",
                               fmt::arg("mg_id", mg.id), fmt::arg("vset_id", mg.vset_id()),
                               fmt::arg("vint_id", mg.vint_id()), fmt::arg("iter_id", iter_id));
        }
        return out;
    };

    std::string gen_code_mg_op(const PlanIR &plan, const VertexSetIR &op) {
        std::optional<MiniGraphIR> mg = plan.get_parent_mg(op);
        if (!mg.has_value()) return gen_code_op(plan, op);
        int dep = op.loop_depth();
        std::optional<VertexSetIR> parent = plan.get_parent_vset(op);
        CHECK(parent.has_value()) << "Logic error (find vset's pruned graph but not its parent)";
        VertexSetIR iter = plan.iter_set.at(op.loop_depth() - 1);
        std::string out;
        if (mg->computed(iter, op)) {
            // Read directly from the pruned graph
            if (op.is_restricted(op.loop_depth())) {
                out += fmt::format("VertexSet s{op_id} = m{mg_id}_adj.bounded(i{dep}_id);\n",
                                   fmt::arg("mg_id", mg->id),
                                   fmt::arg("op_id", op.id),
                                   fmt::arg("iter_id", iter.id),
                                   fmt::arg("dep", dep));
            } else {
                out += fmt::format("VertexSet s{op_id} = m{mg_id}_adj;\n",
                                   fmt::arg("mg_id", mg->id),
                                   fmt::arg("op_id", op.id),
                                   fmt::arg("iter_id", iter.id),
                                   fmt::arg("dep", dep));
            }

            CHECK(!plan.is_last_op(op))
                << "\nLogic error (vset should not be the last op if it can be read directly from a pruned graph)";
        } else if (op.is_edge(op.loop_depth())) {
            std::string upper_bound;
            if (op.is_restricted(op.loop_depth())) {
                upper_bound = fmt::format(", m{mg_id}_adj.vid()", fmt::arg("mg_id", mg->id));
            }
            if (!plan.is_last_op(op)) {
                // intersect and return vertex set
                out += fmt::format(
                        "VertexSet s{op_id} = s{parent_id}.intersect(m{mg_id}_adj{upper_bound});\n",
                        fmt::arg("op_id", op.id),
                        fmt::arg("parent_id", parent->id),
                        fmt::arg("mg_id", mg->id),
                        fmt::arg("upper_bound", upper_bound));
            } else {
                // intersect and return counter
                out += fmt::format("counter += s{parent_id}.intersect_cnt(m{mg_id}_adj{upper_bound});\n",
                                   fmt::arg("parent_id", parent->id),
                                   fmt::arg("mg_id", mg->id),
                                   fmt::arg("upper_bound", upper_bound));
            }
        } else {
            // VertexInduced: subtraction | EdgeInduced: remove one (should not arrive here)
            CHECK(VertexSetIR::adjMatType == AdjMatType::VertexInduced)
                << "Logic error (EdgeInduced patterns does not require set subtraction)";
            std::string upper_bound;
            if (op.is_restricted(op.loop_depth())) {
                upper_bound = fmt::format(", m{mg_id}_adj.vid()", fmt::arg("mg_id", mg->id));
            }
            if (!plan.is_last_op(op)) {
                // subtract and return vertex set
                out += fmt::format(
                        "VertexSet s{op_id} = s{parent_id}.subtract(m{mg_id}_adj{upper_bound});\n",
                        fmt::arg("op_id", op.id),
                        fmt::arg("parent_id", parent->id),
                        fmt::arg("mg_id", mg->id),
                        fmt::arg("upper_bound", upper_bound));
            } else {
                // subtract and return counter
                out += fmt::format("counter += s{parent_id}.subtract_cnt(m{mg_id}_adj{upper_bound});\n",
                                   fmt::arg("parent_id", parent->id),
                                   fmt::arg("mg_id", mg->id),
                                   fmt::arg("upper_bound", upper_bound));
            }
        }
        if (!plan.is_last_op(op))
            out += gen_indent(dep) + fmt::format("if (s{op_id}.size() == 0) continue;\n", fmt::arg("op_id", op.id));
        return out;
    };

    std::string gen_code_iep(const PlanIR &plan, size_t group_id) {
        int val = plan.iep_vals.at(group_id);
        const auto &group = plan.iep_groups.at(group_id);
        std::string out = fmt::format("counter += {}ll", val);
        for (size_t set_id = 0; set_id < group.size(); set_id++) {
            const auto &set = group.at(set_id);
            if (set.size() == 1) {
                int set_id = set.at(0);
                const VertexSetIR &left = plan.iep_set.at(set_id);
                out += fmt::format(" * s{left_id}.size()",
                                   fmt::arg("left_id", left.id));
            } else if (set.size() == 2) {
                const VertexSetIR &left = plan.iep_set.at(set.at(0));
                const VertexSetIR &right = plan.iep_set.at(set.at(1));
                if (left == right) {
                    out += fmt::format(" * s{left_id}.size()",
                                       fmt::arg("left_id", left.id));
                } else {
                    out += fmt::format(" * s{left_id}.intersect_cnt(s{right_id})",
                                       fmt::arg("left_id", left.id),
                                       fmt::arg("right_id", right.id));
                }

            } else {
                const VertexSetIR &left = plan.iep_set.at(set.at(0));
                out += fmt::format(" * s{left_id}", fmt::arg("left_id", left.id));
                for (size_t i = 1; i < set.size() - 1; ++i) {
                    const VertexSetIR &right = plan.iep_set.at(set.at(i));
                    out += fmt::format(".intersect(s{right_id})", fmt::arg("right_id", right.id));
                }
                const VertexSetIR &right = plan.iep_set.at(set.at(set.size() - 1));
                out += fmt::format(".intersect_cnt(s{right_id})", fmt::arg("right_id", right.id));
            }
        }
        out += ";\n";
        return out;
    }

    std::string gen_comment_iep(const PlanIR &plan, size_t group_id) {
        int val = plan.iep_vals.at(group_id);
        std::string group_str, comp_str;
        const auto &group = plan.iep_groups.at(group_id);
        size_t j = 0;
        for (const auto &set: group) {
            size_t i = 0;
            group_str += "(";
            comp_str += "|";
            for (auto set_id: set) {
                group_str += std::to_string(set_id);
                comp_str += fmt::format("VSet({})", plan.iep_set.at(set_id).id);
                if (i++ != set.size() - 1) {
                    group_str += " ";
                    comp_str += " & ";
                }
            }
            group_str += ")";
            comp_str += "|";
            if (j++ != group.size() - 1) {
                group_str += ",";
                comp_str += "*";
            }
        }
        return fmt::format("/* Val: {val} | Group: {group_str} | Comp: {compute_str} */\n",
                           fmt::arg("val", val), fmt::arg("group_str", group_str), fmt::arg("compute_str", comp_str));
    }

    std::string gen_code_omp(PlanIR plan, CodeGenConfig config) {
        Timer t;
        std::ostringstream out;
        if (EnableProfling) out << "#include \"plan_profile.h\"\n";
        else out << "#include \"plan.h\"\n";
        out << "namespace minigraph {\n";
        out << "\tuint64_t pattern_size() {return " << plan.p_size << ";}\n";
        out << "\tvoid plan(const GraphType* graph, Context& ctx){\n";
        if (EnableProfling) out << "\t\tVertexSet::profiler = ctx.profiler;\n";

        switch (config.pruningType) {
            case (PruningType::Eager):
                out << "\t\tusing MiniGraphType = MiniGraphEager;\n";
                break;
            case (PruningType::Static):
                out << "\t\tusing MiniGraphType = MiniGraphLazy;\n";
                break;
            case (PruningType::Online):
                out << "\t\tusing MiniGraphType = MiniGraphOnline;\n";
                break;
            case (PruningType::CostModel):
                out << "\t\tusing MiniGraphType = MiniGraphCostModel;\n";
                break;
            default:
                break;
        }
        if (config.pruningType != PruningType::None) out << "\t\tMiniGraphIF::DATA_GRAPH = graph;\n";
        out << "\t\tVertexSetType::MAX_DEGREE = graph->get_maxdeg();\n";
        out << "#pragma omp parallel num_threads(ctx.num_threads) default(none) shared(ctx, graph)\n\t\t{ // pragma parallel \n";
        out << "\t\t\tcc &counter = ctx.per_thread_result.at(omp_get_thread_num());\n";
        out << "\t\t\tcc &handled = ctx.per_thread_handled.at(omp_get_thread_num());\n";
        out << "\t\t\tdouble start = omp_get_wtime();\n";
        out << "\t\t\tctx.iep_redundency = " << plan.iep_redundancy << ";\n";
        out << "#pragma omp for schedule(dynamic, 1) nowait\n";
        out << "\t\t\tfor (IdType i0_id = 0; i0_id < graph->get_vnum(); i0_id++) { // loop-0 begin\n";
        int max_dep = plan.p_size - 1;
        const auto &set_ops = plan.set_ops;
        switch (config.pruningType) {
            case (PruningType::None):
                if (config.adjMatType != AdjMatType::EdgeInducedIEP || plan.iep_num <= 1) {
                    for (int dep = 0; dep < max_dep; dep++) {
                        // code for reading adj from the graph
                        out << gen_indent(dep) << gen_code_read_adj(plan, dep);
                        // code for computation at this loop
                        const auto &ops = set_ops.at(dep);
                        for (const auto &op: ops) {
                            out << gen_indent(dep) << gen_code_op(plan, op);
                            out << gen_indent(dep) << op;
                        }
                        // code for iterating next loop
                        if (dep == plan.p_size - 2) continue;
                        out << gen_indent(dep) << gen_code_iter(plan, dep);
                    }
                } else {
                    assert(plan.iep_num + plan.iep_depth == plan.p_size - 1);
                    for (int dep = 0; dep < plan.iep_depth; dep++) {
                        // code for reading adj from the graph
                        out << gen_indent(dep) << gen_code_read_adj(plan, dep);
                        // code for computation at this loop
                        const auto &ops = set_ops.at(dep);
                        for (const auto &op: ops) {
                            out << gen_indent(dep) << gen_code_op(plan, op);
                            out << gen_indent(dep) << op;
                        }
                        // code for iterating next loop
                        if (dep == plan.p_size - 2) continue;
                        out << gen_indent(dep) << gen_code_iter(plan, dep);
                    }
                    // Code for IEP
                    int dep = plan.iep_depth;
                    out << gen_indent(dep) << gen_code_read_adj(plan, dep);
                    // code for computation at this loop
                    const auto &ops = set_ops.at(dep);
                    for (const auto &op: ops) {
                        out << gen_indent(dep) << gen_code_op(plan, op);
                        out << gen_indent(dep) << op;
                    }

                    for (size_t group_id = 0; group_id < plan.iep_groups.size(); group_id++) {
                        out << gen_indent(dep) << gen_code_iep(plan, group_id);
                        out << gen_indent(dep) << gen_comment_iep(plan, group_id);
                    }
                }
                break;

            default: // enable pruning
                if (config.adjMatType != AdjMatType::EdgeInducedIEP || plan.iep_num <= 1) {
                    for (int dep = 0; dep < max_dep; dep++) {
                        // code for reading adj from the graph
                        out << gen_indent(dep) << gen_code_read_adj(plan, dep);
                        if (dep > 0) out << gen_code_mg_adj(plan, dep);
                        // code for computation at this loop
                        const auto &ops = set_ops.at(dep);
                        for (const auto &op: ops) {
                            out << gen_indent(dep) << gen_code_mg_op(plan, op);
                            out << gen_indent(dep) << op;
                        }
                        if (dep == plan.p_size - 2) continue;
                        // code for building pruned graphs
                        const auto &mgs = plan.mg_ops.at(dep);
                        for (const auto &mg: mgs) {
                            out << gen_indent(dep) << gen_code_mg_init(plan, mg);
                            out << gen_indent(dep) << mg;
                            out << gen_indent(dep) << gen_code_mg_build(plan, mg);
                        }

                        for (const auto &mg: plan.mg_used.at(dep + 1)) {
                            out << gen_indent(dep) << gen_code_mg_indice(plan, mg, dep);
                        }
                        // code for iterating next loop
                        out << gen_indent(dep) << gen_code_iter(plan, dep);
                    }
                } else {
                    assert(plan.iep_num + plan.iep_depth == plan.p_size - 1);
                    for (int dep = 0; dep < plan.iep_depth; dep++) {
                        // code for reading adj from the graph
                        out << gen_indent(dep) << gen_code_read_adj(plan, dep);
                        if (dep > 0) out << gen_code_mg_adj(plan, dep);
                        // code for computation at this loop
                        const auto &ops = set_ops.at(dep);
                        for (const auto &op: ops) {
                            out << gen_indent(dep) << gen_code_mg_op(plan, op);
                            out << gen_indent(dep) << op;
                        }
                        if (dep == plan.p_size - 2) continue;
                        // code for building pruned graphs
                        const auto &mgs = plan.mg_ops.at(dep);
                        for (const auto &mg: mgs) {
                            out << gen_indent(dep) << gen_code_mg_init(plan, mg);
                            out << gen_indent(dep) << mg;
                            out << gen_indent(dep) << gen_code_mg_build(plan, mg);
                        }

                        for (const auto &mg: plan.mg_used.at(dep + 1)) {
                            out << gen_indent(dep) << gen_code_mg_indice(plan, mg, dep);
                        }
                        // code for iterating next loop
                        out << gen_indent(dep) << gen_code_iter(plan, dep);
                    }
                    int dep = plan.iep_depth;
                    if (dep > 0) out << gen_code_mg_adj(plan, dep);
                    out << gen_indent(dep) << gen_code_read_adj(plan, dep);
                    // code for computation at this loop
                    const auto &ops = set_ops.at(dep);
                    for (const auto &op: ops) {
                        out << gen_indent(dep) << gen_code_mg_op(plan, op);
                        out << gen_indent(dep) << op;
                    }

                    for (size_t group_id = 0; group_id < plan.iep_groups.size(); group_id++) {
                        out << gen_indent(dep) << gen_code_iep(plan, group_id);
                        out << gen_indent(dep) << gen_comment_iep(plan, group_id);
                    }
                }
                break;
        }
        if (plan.iep_num <= 1) {
            for (int dep = max_dep - 1; dep >= 0; dep--) {
                if (dep == 0) out << gen_indent(0) << "handled+=1;\n";
                out << gen_indent(dep) << "} // loop-" << std::to_string(dep) << " end\n";
            }
        } else {
            for (int dep = plan.iep_depth; dep >= 0; dep--) {
                if (dep == 0) out << gen_indent(0) << "handled+=1;\n";
                out << gen_indent(dep) << "} // loop-" << std::to_string(dep) << " end\n";
            }
        };

        out << "\t\t\tctx.per_thread_time.at(omp_get_thread_num()) = omp_get_wtime() - start;\n";
        out << "\t\t} // pragma parallel\n";
        out << "\t} // plan\n";
        out << "} // namespace minigraph \n";

        out << "extern \"C\" void plan(const minigraph::GraphType* graph, minigraph::Context& ctx){return minigraph::plan(graph, ctx);};";
        LOG(INFO) << "Code Generation Time: " << t.Passed() << "s";
        return out.str();
    };

    // all the mingraphs used inside the loop-th loop
    std::vector<MiniGraphIR> gen_used_mg(const PlanIR& plan, const CodeGenConfig& config, int loop){
        std::vector<MiniGraphIR> used_mg;
        if (config.pruningType != PruningType::None && loop > 0) {
            for (int dep = loop; dep < plan.p_size - 1; dep++) {
                for (auto mg: plan.mg_used.at(dep)) {
                    if (mg.loop_depth() < loop) {
                        bool exited = false;
                        for (auto _mg: used_mg) {
                            if (_mg.id == mg.id) exited = true;
                        }
                        if (!exited) used_mg.push_back(mg);
                    }
                }
            }
            for (auto mg: plan.mg_ops.at(loop)) {
                auto parent = plan.get_parent_mg(mg);
                if (parent.has_value() && parent->loop_depth() < loop) {
                    bool exited = false;
                    for (auto _mg: used_mg) {
                        if (_mg.id == parent->id) exited = true;
                    }
                    if (!exited) used_mg.push_back(parent.value());
                }
            }
        }
        return used_mg;
    }

    std::vector<VertexSetIR> gen_used_set(const PlanIR& plan, const CodeGenConfig& config, int loop) {
        std::vector<VertexSetIR> used_set;
        if (loop > 0) {
            const VertexSetIR &iter = plan.iter_set.at(loop - 1);
            for (VertexSetIR set: plan.set_ops.at(loop)) {
                std::optional<VertexSetIR> parent = plan.get_parent_vset(set, loop - 1);
                if (parent.has_value() && parent->id != iter.id) {
                    bool exited = false;
                    for (auto _set: used_set) {
                        if (_set.id == parent->id) exited = true;
                    }
                    if (!exited) used_set.push_back(parent.value());
                }
            }
        }
        return used_set;
    }

    std::set<int> gen_used_adj(const PlanIR& plan, const CodeGenConfig& config, int loop){
        std::set<int> used_adj;
        if (loop > 0) {
            for (int i = loop; i < plan.p_size - 1; ++i) {
                const auto &sets = plan.set_ops.at(i);
                for (auto &set: sets) {
                    if (!plan.get_parent_vset(set, loop - 1).has_value()) {
                        for (int i = 0; i < loop; ++i) {
                            used_adj.insert(i);
                        }
                    }
                }
            }
        }
        return used_adj;
    }

    // loop: the loop at which the next parallel region is evoked
    std::string gen_code_tbb_call(const PlanIR &plan, const CodeGenConfig &config, int loop, int indent_dep) {
        std::ostringstream out;
        if (loop >= plan.get_serial_loop()) return "";
        if (config.parType != ParallelType::Nested && config.parType != ParallelType::NestedRt) return "";
        assert(loop > 0);
        int iter_id = plan.iter_set.at(loop - 1).id;
        const VertexSetIR &iter = plan.iter_set.at(loop - 1);
        std::vector<MiniGraphIR> used_mg = gen_used_mg(plan, config, loop);
        std::vector<VertexSetIR> used_set = gen_used_set(plan, config, loop);
        std::set<int> used_adj = gen_used_adj(plan, config, loop);
        if (config.parType == ParallelType::Nested) {
            out << gen_indent_tbb(indent_dep) + "if (true) ";
        } else if (config.parType == ParallelType::NestedRt) {
            // TODO try using different heuristic logic here
            int factor = 4;
            int avg_deg = plan.meta.num_edge / plan.meta.num_vertex;
            if (plan.meta.max_degree / avg_deg > 100) {
                out << gen_indent_tbb(indent_dep) + fmt::format("if (s{iter_id}.size() > std::min({factor} * {avg_deg}, 100)) ",
                                                                fmt::arg("iter_id", iter_id),
                                                                fmt::arg("avg_deg", avg_deg),
                                                                fmt::arg("factor", factor));
            } else {
                out << gen_indent_tbb(indent_dep) + fmt::format("if (s{iter_id}.size() > {factor} * {avg_deg}) ",
                                                                fmt::arg("iter_id", iter_id),
                                                                fmt::arg("avg_deg", avg_deg),
                                                                fmt::arg("factor", factor));
            }

        }

        int grain_size = 1;

        out << "{\n";
        out << gen_indent_tbb(indent_dep + 1) +
               fmt::format("tbb::parallel_for(tbb::blocked_range<size_t>(0, s{iter_id}.size(), {grain_size}), Loop{dep}",
                           fmt::arg("grain_size", grain_size),
                           fmt::arg("iter_id", iter_id),
                           fmt::arg("dep", loop));
        // Args
        out << "(ctx";
        for (int dep: used_adj) {
            out << fmt::format(", i{}_adj", dep);
        }

        for (auto set: used_set) {
            out << fmt::format(", s{}", set.id);
        }

        out << fmt::format(", s{}", iter_id);

        if (config.pruningType != PruningType::None) {
            for (auto mg: plan.mg_used.at(loop)) {
                if (!skip_build_indices(plan, mg, iter)) out << fmt::format(", m{}_s{}", mg.id, iter_id);
            }

            for (auto mg: used_mg) {
                out << fmt::format(", m{}", mg.id);
            }
        }
        out << "), tbb::auto_partitioner()";
        out << "); continue;\n";
        out << gen_indent_tbb(indent_dep) << "}\n";
        return out.str();
    }

    std::string gen_code_tbb_loop(const PlanIR &plan, const CodeGenConfig &config, int loop) {
        std::ostringstream out;
        int iter_id = -1;
        VertexSetIR iter;
        if (loop > 0) {
            iter_id = plan.iter_set.at(loop - 1).id;
            iter = plan.iter_set.at(loop - 1);
        }
        std::vector<MiniGraphIR> used_mg = gen_used_mg(plan, config, loop);
        std::vector<VertexSetIR> used_set = gen_used_set(plan, config, loop);
        std::set<int> used_adj = gen_used_adj(plan, config, loop);

        out << "\tclass Loop" << loop << "\n\t{\n";

        // Private Variables
        out << "\tprivate:\n";
        out << "\t\tContext& ctx;\n";

        if (!used_adj.empty()) out << "\t\t// Adjacent Lists\n";
        for (int dep: used_adj) {
            out << fmt::format("\t\tVertexSet& i{}_adj;\n", dep);
        }

        if (!used_set.empty()) out << "\t\t// Parent Intermediates\n";
        for (auto set: used_set) {
            out << fmt::format("\t\tVertexSet& s{};\n", set.id);
        }

        if (loop > 0) out << "\t\t// Iterate Set\n" << fmt::format("\t\tVertexSet& s{};\n", iter_id);

        if (config.pruningType != PruningType::None) {
            if (!plan.mg_used.at(loop).empty()) out << "\t\t// MiniGraphs Indices\n";
            for (auto mg: plan.mg_used.at(loop)) {
                if (!skip_build_indices(plan, mg, iter))
                    out << fmt::format("\t\tManagedContainer& m{}_s{};\n", mg.id, iter_id);
            }

            if (!used_mg.empty()) out << "\t\t// MiniGraphs\n";
            for (auto mg: used_mg) {
                std::string mgType = gen_mg_type(plan, mg);
                out << fmt::format("\t\t{}& m{};\n", mgType, mg.id);
            }
        }
        out << "\tpublic:\n";
        // Constructor
        // Args
        out << "\t\tLoop" << loop << "(Context& _ctx";

        for (int dep: used_adj) {
            out << fmt::format(", VertexSet& _i{}_adj", dep);
        }

        for (auto set: used_set) {
            out << fmt::format(", VertexSet& _s{}", set.id);
        }

        if (loop > 0) out << ", VertexSet& _s" << iter_id;

        if (config.pruningType != PruningType::None) {
            for (auto mg: plan.mg_used.at(loop)) {
                if (!skip_build_indices(plan, mg, iter))
                    out << fmt::format(", ManagedContainer& _m{}_s{}", mg.id, iter_id);
            }

            for (auto mg: used_mg) {
                std::string mgType = gen_mg_type(plan, mg);
                out << fmt::format(", {}& _m{}", mgType, mg.id);
            }
        }
        out << ")";

        // Initialization
        out << ":ctx{_ctx}";

        for (int dep: used_adj) {
            out << fmt::format(", i{}_adj", dep) << "{" << fmt::format("_i{}_adj", dep) << "}";
        }

        for (auto set: used_set) {
            out << fmt::format(", s{}", set.id) << "{" << fmt::format("_s{}", set.id) << "}";
        }

        if (loop > 0) out << fmt::format(", s{}", iter_id) << "{" << fmt::format("_s{}", iter_id) << "}";

        if (config.pruningType != PruningType::None) {
            for (auto mg: plan.mg_used.at(loop)) {
                if (!skip_build_indices(plan, mg, iter))
                    out << fmt::format(", m{}_s{}", mg.id, iter_id) << "{" << fmt::format(" _m{}_s{}", mg.id, iter_id)
                        << "}";
            }

            for (auto mg: used_mg) {
                out << fmt::format(", m{}", mg.id) << "{" << fmt::format("_m{}", mg.id) << "}";
            }
        }
        out << " {};\n";

        // Operator
        out << "\t\tvoid operator()(const tbb::blocked_range<size_t> &r) const {// operator begin\n";
        out << "\t\t\tconst int worker_id = tbb::this_task_arena::current_thread_index();\n";
        out << "\t\t\tcc& counter = ctx.per_thread_result.at(worker_id);\n";
        if (loop > 0) {
            out << "\t\t\t" << fmt::format("for (size_t i{loop}_idx = r.begin(); i{loop}_idx < r.end(); i{loop}_idx++)",
                                           fmt::arg("loop", loop));
        } else {
//            out << "\t\t\tcc& handled = ctx.per_thread_handled.at(worker_id);\n";
//            out << "\t\t\t" << "double& time = ctx.per_thread_time.at(worker_id);\n";
//            out << "\t\t\t" << "tick_count t1 = tick_count::now();\n";
            out << "\t\t\t" << fmt::format("for (size_t i{loop}_id = r.begin(); i{loop}_id < r.end(); i{loop}_id++)",
                                           fmt::arg("loop", loop));
        }
        out << " { // loop-" << loop << "begin\n";
        int max_dep = plan.p_size - 1;
        const auto &set_ops = plan.set_ops;
        switch (config.pruningType) {
            case (PruningType::None):
                if (config.adjMatType != AdjMatType::EdgeInducedIEP || plan.iep_num <= 1) {
                    for (int dep = loop; dep < max_dep; dep++) {
                        int indent_dep = dep - loop;
                        // code for reading adj from the graph
                        out << gen_indent_tbb(indent_dep) << gen_code_read_adj(plan, dep);
                        // code for computation at this loop
                        const auto &ops = set_ops.at(dep);
                        for (const auto &op: ops) {
                            out << gen_indent_tbb(indent_dep) << gen_code_op(plan, op);
                            out << gen_indent_tbb(indent_dep) << op;
                        }
                        // skip iterating next loop
                        if (dep == plan.p_size - 2) continue;

                        // code for calling parallel nested loop
                        out << gen_code_tbb_call(plan, config, dep + 1, indent_dep);

                        // code for serial executing next loop
                        out << gen_indent_tbb(indent_dep) << gen_code_iter(plan, dep);
                    }
                } else {
                    assert(plan.iep_num + plan.iep_depth == plan.p_size - 1);
                    for (int dep = loop; dep < plan.iep_depth; dep++) {
                        int indent_dep = dep - loop;

                        // code for reading adj from the graph
                        out << gen_indent_tbb(indent_dep) << gen_code_read_adj(plan, dep);
                        // code for computation at this loop
                        const auto &ops = set_ops.at(dep);
                        for (const auto &op: ops) {
                            out << gen_indent_tbb(indent_dep) << gen_code_op(plan, op);
                            out << gen_indent_tbb(indent_dep) << op;
                        }

                        // code for iterating next loop
                        if (dep == plan.p_size - 2) continue;

                        // code for calling parallel nested loop
                        out << gen_code_tbb_call(plan, config, dep + 1, indent_dep);

                        out << gen_indent_tbb(indent_dep) << gen_code_iter(plan, dep);
                    }
                    // Code for IEP
                    int dep = plan.iep_depth;
                    int indent_dep = dep - loop;

                    out << gen_indent_tbb(indent_dep) << gen_code_read_adj(plan, dep);
                    // code for computation at this loop
                    const auto &ops = set_ops.at(dep);
                    for (const auto &op: ops) {
                        out << gen_indent_tbb(indent_dep) << gen_code_op(plan, op);
                        out << gen_indent_tbb(indent_dep) << op;
                    }

                    for (size_t group_id = 0; group_id < plan.iep_groups.size(); group_id++) {
                        out << gen_indent_tbb(indent_dep) << gen_code_iep(plan, group_id);
                        out << gen_indent_tbb(indent_dep) << gen_comment_iep(plan, group_id);
                    }
                }
                break;

            default: // enable pruning
                if (config.adjMatType != AdjMatType::EdgeInducedIEP || plan.iep_num <= 1) {
                    for (int dep = loop; dep < max_dep; dep++) {
                        // code for reading adj from the graph
                        int indent_dep = dep - loop;
                        out << gen_indent_tbb(indent_dep) << gen_code_read_adj(plan, dep);
                        if (dep > 0) out << gen_code_mg_adj(plan, dep, indent_dep);
                        // code for computation at this loop
                        const auto &ops = set_ops.at(dep);
                        for (const auto &op: ops) {
                            out << gen_indent_tbb(indent_dep) << gen_code_mg_op(plan, op);
                            out << gen_indent_tbb(indent_dep) << op;
                        }
                        if (dep == plan.p_size - 2) continue;
                        // code for building pruned graphs
                        const auto &mgs = plan.mg_ops.at(dep);
                        for (const auto &mg: mgs) {
                            out << gen_indent_tbb(indent_dep) << gen_code_mg_init(plan, mg);
                            out << gen_indent_tbb(indent_dep) << mg;
                            out << gen_indent_tbb(indent_dep) << gen_code_mg_build(plan, mg);
                        }

                        for (const auto &mg: plan.mg_used.at(dep + 1)) {
                            out << gen_indent_tbb(indent_dep) << gen_code_mg_indice(plan, mg, dep);
                        }

                        // code for calling parallel nested loop
                        out << gen_code_tbb_call(plan, config, dep + 1, indent_dep);

                        // code for serially iterating next loop
                        out << gen_indent_tbb(indent_dep) << gen_code_iter(plan, dep);
                    }
                } else {
                    assert(plan.iep_num + plan.iep_depth == plan.p_size - 1);
                    for (int dep = loop; dep < plan.iep_depth; dep++) {
                        int indent_dep = dep - loop;
                        // code for reading adj from the graph
                        out << gen_indent_tbb(indent_dep) << gen_code_read_adj(plan, dep);
                        if (dep > 0) out << gen_code_mg_adj(plan, dep, indent_dep);
                        // code for computation at this loop
                        const auto &ops = set_ops.at(dep);
                        for (const auto &op: ops) {
                            out << gen_indent_tbb(indent_dep) << gen_code_mg_op(plan, op);
                            out << gen_indent_tbb(indent_dep) << op;
                        }
                        if (dep == plan.p_size - 2) continue;
                        // code for building pruned graphs
                        const auto &mgs = plan.mg_ops.at(dep);
                        for (const auto &mg: mgs) {
                            out << gen_indent_tbb(indent_dep) << gen_code_mg_init(plan, mg);
                            out << gen_indent_tbb(indent_dep) << mg;
                            out << gen_indent_tbb(indent_dep) << gen_code_mg_build(plan, mg);
                        }

                        for (const auto &mg: plan.mg_used.at(dep + 1)) {
                            out << gen_indent_tbb(indent_dep) << gen_code_mg_indice(plan, mg, dep);
                        }

                        // code for calling parallel nested loop
                        out << gen_code_tbb_call(plan, config, dep + 1, indent_dep);

                        // code for iterating next loop
                        out << gen_indent_tbb(indent_dep) << gen_code_iter(plan, dep);
                    }
                    int dep = plan.iep_depth;
                    int indent_dep = dep - loop;
                    if (dep > 0) out << gen_code_mg_adj(plan, dep, indent_dep);
                    out << gen_indent_tbb(indent_dep) << gen_code_read_adj(plan, dep);
                    // code for computation at this loop
                    const auto &ops = set_ops.at(dep);
                    for (const auto &op: ops) {
                        out << gen_indent_tbb(indent_dep) << gen_code_mg_op(plan, op);
                        out << gen_indent_tbb(indent_dep) << op;
                    }

                    for (size_t group_id = 0; group_id < plan.iep_groups.size(); group_id++) {
                        out << gen_indent_tbb(indent_dep) << gen_code_iep(plan, group_id);
                        out << gen_indent_tbb(indent_dep) << gen_comment_iep(plan, group_id);
                    }
                }
                break;
        }
        if (plan.iep_num <= 1) {
            for (int dep = max_dep - 1; dep >= loop; dep--) {
                int indent_dep = dep - loop;
//                if (dep == 0 && loop == 0) out << gen_indent(dep) << "handled += 1;\n";
                out << gen_indent_tbb(indent_dep) << "} // loop-" << std::to_string(dep) << " end\n";
            }
        } else {
            for (int dep = plan.iep_depth; dep >= loop; dep--) {
                int indent_dep = dep - loop;
//                if (dep == 0 && loop == 0) out << gen_indent(dep) << "handled += 1;\n";
                out << gen_indent_tbb(indent_dep) << "} // loop-" << std::to_string(dep) << " end\n";
            }
        };

        // if (loop <= 2) out << "\t\t\tctx.per_thread_tick.at(worker_id) = tick_count::now();\n";
        out << "\t\t} // operator end\n";
        out << "\t}; // Loop\n\n";
        return out.str();
    }

    std::string gen_code_nested(PlanIR plan, CodeGenConfig config) {
        Timer t;
        std::ostringstream out;
        if (EnableProfling) out << "#include \"plan_profile.h\"\n";
        else out << "#include \"plan.h\"\n";
        // out << "#include \"oneapi/tbb/parallel_for.h\"\n";
        out << "namespace minigraph {\n";
        out << "\tuint64_t pattern_size() {return " << plan.p_size << ";}\n";
        out << "\tstatic const Graph * graph;\n";

        switch (config.pruningType) {
            case (PruningType::Eager):
                out << "\tusing MiniGraphType = MiniGraphEager;\n";
                break;
            case (PruningType::Static):
                out << "\tusing MiniGraphType = MiniGraphLazy;\n";
                break;
            case (PruningType::Online):
                out << "\tusing MiniGraphType = MiniGraphOnline;\n";
                break;
            case (PruningType::CostModel):
                out << "\tusing MiniGraphType = MiniGraphCostModel;\n";
                break;
            default:
                break;
        }
        for (int loop = plan.get_serial_loop() - 1; loop >= 0; loop--) {
            out << gen_code_tbb_loop(plan, config, loop);
        }
        out << "\tvoid plan(const GraphType* _graph, Context& ctx){ // plan \n";
        if (EnableProfling) {
            out << "\t\tVertexSet::profiler = ctx.profiler;\n";
        }
        out << "\t\tctx.tick_begin = tbb::tick_count::now();\n";
        out << "\t\tctx.iep_redundency = " << plan.iep_redundancy << ";\n";
        out << "\t\tgraph = _graph;\n";
        if (config.pruningType != PruningType::None) out << "\t\tMiniGraphIF::DATA_GRAPH = graph;\n";
        out << "\t\tVertexSetType::MAX_DEGREE = graph->get_maxdeg();\n";
        out << "\t\ttbb::parallel_for(tbb::blocked_range<size_t>(0, graph->get_vnum()), Loop0(ctx), tbb::simple_partitioner());\n";
        out << "\t} // plan\n";
        out << "} // minigraph\n";
        return out.str();
    }

    std::string gen_code(const std::string &adj_mat, CodeGenConfig config, MetaData meta) {
        VertexSetIR::adjMatType = config.adjMatType;
        PlanIR plan = create_plan(adj_mat, config, meta);
        CurConfig = config;
        if (config.pruningType != PruningType::None) {
            plan = create_plan_mg(plan, config);
        }
        if (config.runnerType == RunnerType::Profiling) {
            EnableProfling = true;
        } else {
            EnableProfling = false;
        }
        if (config.parType == ParallelType::OpenMP) {
            return gen_code_omp(plan, config);
        } else if (config.parType == ParallelType::TbbTop
                   || config.parType == ParallelType::Nested
                   || config.parType == ParallelType::NestedRt) {
            return gen_code_nested(plan, config);
        } else {
            return ""; // TODO distributed
        }
    };
}