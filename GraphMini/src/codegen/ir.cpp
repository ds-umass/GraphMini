//
// Created by ubuntu on 1/1/23.
//

#include "ir.h"
#include <iostream>

namespace minigraph
{

    VertexSetIR::VertexSetIR(const EdgeIR &edges, const EdgeRestrictIR restricts, int loop_depth) :
            m_loop_depth{loop_depth} {
        for (int i = 0; i <= loop_depth; i++) {
            m_edges[i] = edges[i];
            m_restricts[i] = restricts[i];
        };
    }

    bool VertexSetIR::share_at_least_one_parent_node(const VertexSetIR &other) const {
        int max_dep = std::min(other.loop_depth(), loop_depth());
        for (int dep = 0; dep <= max_dep; dep++){
            if (is_edge(dep) && other.is_edge(dep)) return true;
        }
        return false;
    }

    bool VertexSetIR::operator==(const VertexSetIR &rhs) const {
        return m_edges == rhs.m_edges &&
               m_restricts == rhs.m_restricts &&
               m_loop_depth == rhs.m_loop_depth;
    }

    bool VertexSetIR::same_iep_computation(const VertexSetIR &rhs) const {
        switch (adjMatType) {
            case AdjMatType::VertexInduced:
                return false;

            default:
                int min_dep = std::min(rhs.loop_depth(), loop_depth());
                for (int i = 0; i <= min_dep; ++i) {
                    if (is_edge(i) != rhs.is_edge(i) || is_restricted(i) != rhs.is_restricted(i)) return false;
                }
                return true;
        }
    }

    bool VertexSetIR::operator<(const VertexSetIR &rhs) const {
        if (m_loop_depth < rhs.loop_depth()) return true;
        else if (edge_num() < rhs.edge_num()) return true;
        else if (restrict_num() < rhs.restrict_num()) return true;
        for (int i = 0; i <= m_loop_depth; i++) {
            if (is_edge(i) && !rhs.is_edge(i)) {
                return false;
            }
            if (!is_edge(i) && rhs.is_edge(i)) {
                return true;
            }
        }
        for (int depth = 0; depth <= m_loop_depth; depth++) {
            if (is_restricted(depth) && !rhs.is_restricted(depth)) {
                return false;
            }
            if (!is_restricted(depth) && rhs.is_restricted(depth)) {
                return true;
            }
        }
        return false;
    }

    bool VertexSetIR::is_superset_of(const VertexSetIR &rhs) const {
        if (m_loop_depth > rhs.m_loop_depth) return false;
        for (int i = 0; i <= m_loop_depth; i++) {
            if (adjMatType == AdjMatType::VertexInduced) {
                if ((rhs.is_edge(i) != is_edge(i)) ||
                    (!rhs.is_restricted(i) && is_restricted(i))) {
                    // if any in_edge relation doesn't match or
                    // any restrict in "this" is not presented in "rhs"
                    return false;
                }
            } else {
                if ((!rhs.is_edge(i) && is_edge(i)) ||
                    (!rhs.is_restricted(i) && is_restricted(i))) {
                    // if any in_edge relation doesn't match or
                    // any restrict in "this" is not presented in "rhs"
                    return false;
                }
            }
        }
        return true;
    }

    std::ostream &operator<<(std::ostream &out, const VertexSetIR &v){
        out << "/* VSet(" << v.id << ", " << v.loop_depth() << ") In-Edges: ";
        for (int dep = 0; dep <= v.loop_depth(); dep++){
            if (v.m_edges[dep]) out << dep << " ";
        }
        out << "Restricts: ";
        for (int dep = 0; dep <= v.loop_depth(); dep++){
            if (v.m_restricts[dep]) out << dep << " ";
        }
        out << "*/\n";
        return out;
    };

    bool MiniGraphIR::operator==(const MiniGraphIR &rhs) const {
        return m_intersect == rhs.m_intersect &&
               m_vertices == rhs.m_vertices;
    }

    bool MiniGraphIR::operator<(const MiniGraphIR &rhs) const {
        if (rhs == *this) return false;
        else if (m_vertices < rhs.m_vertices) return true;
        else if (rhs.m_vertices < m_vertices) return false;
        else if (m_intersect < rhs.m_intersect) return true;
        else return false;
    }

    bool MiniGraphIR::is_superset_of(const MiniGraphIR &rhs) const {
        return m_vertices.is_superset_of(rhs.m_vertices) &&
               m_intersect.is_superset_of(rhs.m_intersect);
    }

    bool MiniGraphIR::is_superset_of(const VertexSetIR &vertices, const VertexSetIR &intersect) const {
        return m_vertices.is_superset_of(vertices) &&
               m_intersect.is_superset_of(intersect);
    }

    bool MiniGraphIR::computed(const VertexSetIR &vertices, const VertexSetIR &intersect) const {
        bool flag = is_superset_of(vertices, intersect) &&
                    intersect.loop_depth() - m_intersect.loop_depth() == 1 &&
                    intersect.edge_num() - m_intersect.edge_num() == 1;

        if (!flag) return false;
        for (int i = 0; i <= m_intersect.loop_depth(); ++i) {
            if (m_intersect.is_edge(i) != intersect.is_edge(i)) return false;
            if (m_intersect.is_restricted(i) != intersect.is_restricted(i)) return false;
        }
        return true;
    }

    std::optional<VertexSetIR> PlanIR::get_parent_vset(const VertexSetIR& vset) const {
        std::optional<VertexSetIR> out;
        // vs is op bounded by id
        for (const auto &op: set_ops.at(vset.loop_depth())) {
            if (op.is_superset_of(vset) &&
                op.id != vset.id &&
                op.edge_num() == vset.edge_num() &&
                vset.restrict_num() - op.restrict_num() == 1 &&
                vset.is_restricted(vset.loop_depth()) &&
                !op.is_restricted(vset.loop_depth())) {
                out = op;
            }
            if (out.has_value()) return out;
        }
        // op is vs parent prefix
        if (vset.loop_depth() > 0) { // looking for prefix vertex
            const auto &ops = set_ops.at(vset.loop_depth() - 1);
            for (const auto &op: ops) {
                if (op.is_superset_of(vset)
                    && op.id != vset.id && (!out.has_value() || out->is_superset_of(op)))
                    out = op;
            }
        }
        return out;
    }

    std::optional<VertexSetIR> PlanIR::get_parent_vset(const VertexSetIR& vset, size_t dep) const {
        std::optional<VertexSetIR> out;
        // op is vset's parent prefix
        if (vset.loop_depth() > 0) { // looking for prefix vertex
            const auto &ops = set_ops.at(dep);
            for (const auto &op: ops) {
                if (op.is_superset_of(vset) && op.id != vset.id
                    && (!out.has_value() || out->is_superset_of(op)))
                    out = op;
            }
        }
        return out;
    }
    std::optional<MiniGraphIR> PlanIR::get_parent_mg(const VertexSetIR& intersect) const {
        std::optional<MiniGraphIR> out;
        if(intersect.loop_depth() == 0) return out;
        if(VertexSetIR::adjMatType != AdjMatType::VertexInduced
           && !intersect.is_edge(intersect.loop_depth())) return out;

        VertexSetIR vertices = iter_set.at(intersect.loop_depth() - 1);
        for (int dep = 0; dep < intersect.loop_depth(); ++dep) {
            for (const auto& mg_op : mg_ops.at(dep)) {
                if (!mg_op.is_superset_of(vertices, intersect)) continue;
                if (out.has_value() && out->is_superset_of(mg_op)){
                    out = mg_op;
                } else {
                    out = mg_op;
                }
            }
        }
        return out;
    };
    std::optional<MiniGraphIR> PlanIR::get_parent_mg(const MiniGraphIR& child) const {
        std::optional<MiniGraphIR> out;
        for (int dep = 0; dep <= child.loop_depth(); ++dep) {
            for (const auto& parent : mg_ops.at(dep)) {
                if(!parent.is_superset_of(child) || parent.id == child.id) continue;
                if(is_bounded(parent) && !is_bounded(child)) continue;
                if (out.has_value() && out->is_superset_of(parent)) {
                    out = parent;
                } else {
                    out = parent;
                }
            }
        }
        return out;
    };
    bool PlanIR::is_bounded(const MiniGraphIR &mg) const {
        return mg_bounded.at(mg.id);
    };
    bool PlanIR::is_bounded(const VertexSetIR &vset) const {
        return vset.is_restricted(vset.loop_depth());
    };
    bool PlanIR::is_last_op(const VertexSetIR &op) const {
        return op.loop_depth() == p_size - 2;
    };

    bool PlanIR::is_par(const MiniGraphIR& mg) const {
        return false;
    }
    std::ostream &operator<<(std::ostream &out, const MiniGraphIR &mg){
        VertexSetIR v = mg.m_vertices;
        out << "/* Vertices = VSet(" << v.id << ") In-Edges: ";
        for (int dep = 0; dep <= v.loop_depth(); dep++){
            if (v.m_edges[dep]) out << dep << " ";
        }
        out << "Restricts: ";
        for (int dep = 0; dep <= v.loop_depth(); dep++){
            if (v.m_restricts[dep]) out << dep << " ";
        }
        out << " | Intersect = ";
        v = mg.m_intersect;
        out << "VSet(" << v.id << ") In-Edges: ";
        for (int dep = 0; dep <= v.loop_depth(); dep++){
            if (v.m_edges[dep]) out << dep << " ";
        }
        out << "Restricts: ";
        for (int dep = 0; dep <= v.loop_depth(); dep++){
            if (v.m_restricts[dep]) out << dep << " ";
        }
        out << "*/\n";
        return out;
    };
}