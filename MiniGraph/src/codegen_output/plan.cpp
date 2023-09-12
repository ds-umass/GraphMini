#include "plan.h"
namespace minigraph {
	uint64_t pattern_size() {return 7;}
	static const Graph * graph;
	using MiniGraphType = MiniGraphCostModel;
	class Loop3
	{
	private:
		Context& ctx;
		// Parent Intermediates
		VertexSet& s4;
		// Iterate Set
		VertexSet& s5;
		// MiniGraphs Indices
		// MiniGraphs
		MiniGraphEager& m4;
	public:
		Loop3(Context& _ctx, VertexSet& _s4, VertexSet& _s5, MiniGraphEager& _m4):ctx{_ctx}, s4{_s4}, s5{_s5}, m4{_m4} {};
		void operator()(const tbb::blocked_range<size_t> &r) const {// operator begin
			const int worker_id = tbb::this_task_arena::current_thread_index();
			cc& counter = ctx.per_thread_result.at(worker_id);
			for (size_t i3_idx = r.begin(); i3_idx < r.end(); i3_idx++) { // loop-3begin
				const IdType i3_id = s5[i3_idx];
							VertexSet m4_adj = m4.N(i3_idx);
				VertexSet s6 = m4_adj;
							if (s6.size() == 0) continue;
				/* VSet(6, 3) In-Edges: 0 1 2 3 Restricts: */
				VertexSet s7 = s5.intersect(m4_adj, m4_adj.vid());
							if (s7.size() == 0) continue;
				/* VSet(7, 3) In-Edges: 0 1 2 3 Restricts: 0 1 2 3 */
				auto m4_s7 = m4.indices(s7);
				for (size_t i4_idx = 0; i4_idx < s7.size(); i4_idx++) { // loop-4 begin
								VertexSet m4_adj = m4.N(m4_s7[i4_idx]);
					const IdType i4_id = s7[i4_idx];
					VertexSet s8 = s6.intersect(m4_adj);
								if (s8.size() == 0) continue;
					/* VSet(8, 4) In-Edges: 0 1 2 3 4 Restricts: */
					counter += 1ll * s8.size() * s8.size();
					/* Val: 1 | Group: (0),(1) | Comp: |VSet(8)|*|VSet(8)| */
					counter += -1ll * s8.size();
					/* Val: -1 | Group: (0 1) | Comp: |VSet(8) & VSet(8)| */
					} // loop-4 end
				} // loop-3 end
		} // operator end
	}; // Loop

	class Loop2
	{
	private:
		Context& ctx;
		// Parent Intermediates
		VertexSet& s2;
		// Iterate Set
		VertexSet& s3;
		// MiniGraphs Indices
		// MiniGraphs
		MiniGraphEager& m2;
		MiniGraphEager& m3;
	public:
		Loop2(Context& _ctx, VertexSet& _s2, VertexSet& _s3, MiniGraphEager& _m2, MiniGraphEager& _m3):ctx{_ctx}, s2{_s2}, s3{_s3}, m2{_m2}, m3{_m3} {};
		void operator()(const tbb::blocked_range<size_t> &r) const {// operator begin
			const int worker_id = tbb::this_task_arena::current_thread_index();
			cc& counter = ctx.per_thread_result.at(worker_id);
			for (size_t i2_idx = r.begin(); i2_idx < r.end(); i2_idx++) { // loop-2begin
				const IdType i2_id = s3[i2_idx];
						VertexSet m2_adj = m2.N(i2_idx);
						VertexSet m3_adj = m3.N(i2_idx);
				VertexSet s4 = m2_adj;
						if (s4.size() == 0) continue;
				/* VSet(4, 2) In-Edges: 0 1 2 Restricts: */
				VertexSet s5 = m3_adj.bounded(i2_id);
						if (s5.size() == 0) continue;
				/* VSet(5, 2) In-Edges: 0 1 2 Restricts: 0 1 2 */
				MiniGraphEager m4(false,false);
				/* Vertices = VSet(5) In-Edges: 0 1 2 Restricts: 0 1 2  | Intersect = VSet(4) In-Edges: 0 1 2 Restricts: */
				m4.build(&m2, s5, s4, s5);
				//skip building indices for m4 because they can be obtained directly
				if (s5.size() > 4 * 28) {
					tbb::parallel_for(tbb::blocked_range<size_t>(0, s5.size(), 1), Loop3(ctx, s4, s5, m4), tbb::auto_partitioner()); continue;
				}
				for (size_t i3_idx = 0; i3_idx < s5.size(); i3_idx++) { // loop-3 begin
					const IdType i3_id = s5[i3_idx];
							VertexSet m4_adj = m4.N(i3_idx);
					VertexSet s6 = m4_adj;
							if (s6.size() == 0) continue;
					/* VSet(6, 3) In-Edges: 0 1 2 3 Restricts: */
					VertexSet s7 = s5.intersect(m4_adj, m4_adj.vid());
							if (s7.size() == 0) continue;
					/* VSet(7, 3) In-Edges: 0 1 2 3 Restricts: 0 1 2 3 */
					auto m4_s7 = m4.indices(s7);
					for (size_t i4_idx = 0; i4_idx < s7.size(); i4_idx++) { // loop-4 begin
								VertexSet m4_adj = m4.N(m4_s7[i4_idx]);
						const IdType i4_id = s7[i4_idx];
						VertexSet s8 = s6.intersect(m4_adj);
								if (s8.size() == 0) continue;
						/* VSet(8, 4) In-Edges: 0 1 2 3 4 Restricts: */
						counter += 1ll * s8.size() * s8.size();
						/* Val: 1 | Group: (0),(1) | Comp: |VSet(8)|*|VSet(8)| */
						counter += -1ll * s8.size();
						/* Val: -1 | Group: (0 1) | Comp: |VSet(8) & VSet(8)| */
						} // loop-4 end
					} // loop-3 end
				} // loop-2 end
		} // operator end
	}; // Loop

	class Loop1
	{
	private:
		Context& ctx;
		// Parent Intermediates
		VertexSet& s0;
		// Iterate Set
		VertexSet& s1;
		// MiniGraphs Indices
		// MiniGraphs
		MiniGraphEager& m0;
		MiniGraphEager& m1;
	public:
		Loop1(Context& _ctx, VertexSet& _s0, VertexSet& _s1, MiniGraphEager& _m0, MiniGraphEager& _m1):ctx{_ctx}, s0{_s0}, s1{_s1}, m0{_m0}, m1{_m1} {};
		void operator()(const tbb::blocked_range<size_t> &r) const {// operator begin
			const int worker_id = tbb::this_task_arena::current_thread_index();
			cc& counter = ctx.per_thread_result.at(worker_id);
			for (size_t i1_idx = r.begin(); i1_idx < r.end(); i1_idx++) { // loop-1begin
				const IdType i1_id = s1[i1_idx];
					VertexSet m0_adj = m0.N(i1_idx);
					VertexSet m1_adj = m1.N(i1_idx);
				VertexSet s2 = m0_adj;
					if (s2.size() == 0) continue;
				/* VSet(2, 1) In-Edges: 0 1 Restricts: */
				VertexSet s3 = m1_adj.bounded(i1_id);
					if (s3.size() == 0) continue;
				/* VSet(3, 1) In-Edges: 0 1 Restricts: 0 1 */
				MiniGraphEager m2(false,false);
				/* Vertices = VSet(3) In-Edges: 0 1 Restricts: 0 1  | Intersect = VSet(2) In-Edges: 0 1 Restricts: */
				m2.build(&m0, s3, s2, s3);
				MiniGraphEager m3(true,false);
				/* Vertices = VSet(3) In-Edges: 0 1 Restricts: 0 1  | Intersect = VSet(3) In-Edges: 0 1 Restricts: 0 1 */
				m3.build(&m2, s3, s3, s3);
				//skip building indices for m2 because they can be obtained directly
				//skip building indices for m3 because they can be obtained directly
				if (s3.size() > 4 * 28) {
					tbb::parallel_for(tbb::blocked_range<size_t>(0, s3.size(), 1), Loop2(ctx, s2, s3, m2, m3), tbb::auto_partitioner()); continue;
				}
				for (size_t i2_idx = 0; i2_idx < s3.size(); i2_idx++) { // loop-2 begin
					const IdType i2_id = s3[i2_idx];
						VertexSet m2_adj = m2.N(i2_idx);
						VertexSet m3_adj = m3.N(i2_idx);
					VertexSet s4 = m2_adj;
						if (s4.size() == 0) continue;
					/* VSet(4, 2) In-Edges: 0 1 2 Restricts: */
					VertexSet s5 = m3_adj.bounded(i2_id);
						if (s5.size() == 0) continue;
					/* VSet(5, 2) In-Edges: 0 1 2 Restricts: 0 1 2 */
					MiniGraphEager m4(false,false);
					/* Vertices = VSet(5) In-Edges: 0 1 2 Restricts: 0 1 2  | Intersect = VSet(4) In-Edges: 0 1 2 Restricts: */
					m4.build(&m2, s5, s4, s5);
					//skip building indices for m4 because they can be obtained directly
					if (s5.size() > 4 * 28) {
						tbb::parallel_for(tbb::blocked_range<size_t>(0, s5.size(), 1), Loop3(ctx, s4, s5, m4), tbb::auto_partitioner()); continue;
					}
					for (size_t i3_idx = 0; i3_idx < s5.size(); i3_idx++) { // loop-3 begin
						const IdType i3_id = s5[i3_idx];
							VertexSet m4_adj = m4.N(i3_idx);
						VertexSet s6 = m4_adj;
							if (s6.size() == 0) continue;
						/* VSet(6, 3) In-Edges: 0 1 2 3 Restricts: */
						VertexSet s7 = s5.intersect(m4_adj, m4_adj.vid());
							if (s7.size() == 0) continue;
						/* VSet(7, 3) In-Edges: 0 1 2 3 Restricts: 0 1 2 3 */
						auto m4_s7 = m4.indices(s7);
						for (size_t i4_idx = 0; i4_idx < s7.size(); i4_idx++) { // loop-4 begin
								VertexSet m4_adj = m4.N(m4_s7[i4_idx]);
							const IdType i4_id = s7[i4_idx];
							VertexSet s8 = s6.intersect(m4_adj);
								if (s8.size() == 0) continue;
							/* VSet(8, 4) In-Edges: 0 1 2 3 4 Restricts: */
							counter += 1ll * s8.size() * s8.size();
							/* Val: 1 | Group: (0),(1) | Comp: |VSet(8)|*|VSet(8)| */
							counter += -1ll * s8.size();
							/* Val: -1 | Group: (0 1) | Comp: |VSet(8) & VSet(8)| */
							} // loop-4 end
						} // loop-3 end
					} // loop-2 end
				} // loop-1 end
		} // operator end
	}; // Loop

	class Loop0
	{
	private:
		Context& ctx;
	public:
		Loop0(Context& _ctx):ctx{_ctx} {};
		void operator()(const tbb::blocked_range<size_t> &r) const {// operator begin
			const int worker_id = tbb::this_task_arena::current_thread_index();
			cc& counter = ctx.per_thread_result.at(worker_id);
			for (size_t i0_id = r.begin(); i0_id < r.end(); i0_id++) { // loop-0begin
				VertexSet i0_adj = graph->N(i0_id);
				VertexSet s0 = i0_adj;
				if (s0.size() == 0) continue;
				/* VSet(0, 0) In-Edges: 0 Restricts: */
				VertexSet s1 = s0.bounded(i0_id);
				/* VSet(1, 0) In-Edges: 0 Restricts: 0 */
				MiniGraphEager m0(false,false);
				/* Vertices = VSet(1) In-Edges: 0 Restricts: 0  | Intersect = VSet(0) In-Edges: 0 Restricts: */
				m0.build(s1, s0, s1);
				MiniGraphEager m1(true,false);
				/* Vertices = VSet(1) In-Edges: 0 Restricts: 0  | Intersect = VSet(1) In-Edges: 0 Restricts: 0 */
				m1.build(&m0, s1, s1, s1);
				//skip building indices for m0 because they can be obtained directly
				//skip building indices for m1 because they can be obtained directly
				if (s1.size() > 4 * 28) {
					tbb::parallel_for(tbb::blocked_range<size_t>(0, s1.size(), 1), Loop1(ctx, s0, s1, m0, m1), tbb::auto_partitioner()); continue;
				}
				for (size_t i1_idx = 0; i1_idx < s1.size(); i1_idx++) { // loop-1 begin
					const IdType i1_id = s1[i1_idx];
					VertexSet m0_adj = m0.N(i1_idx);
					VertexSet m1_adj = m1.N(i1_idx);
					VertexSet s2 = m0_adj;
					if (s2.size() == 0) continue;
					/* VSet(2, 1) In-Edges: 0 1 Restricts: */
					VertexSet s3 = m1_adj.bounded(i1_id);
					if (s3.size() == 0) continue;
					/* VSet(3, 1) In-Edges: 0 1 Restricts: 0 1 */
					MiniGraphEager m2(false,false);
					/* Vertices = VSet(3) In-Edges: 0 1 Restricts: 0 1  | Intersect = VSet(2) In-Edges: 0 1 Restricts: */
					m2.build(&m0, s3, s2, s3);
					MiniGraphEager m3(true,false);
					/* Vertices = VSet(3) In-Edges: 0 1 Restricts: 0 1  | Intersect = VSet(3) In-Edges: 0 1 Restricts: 0 1 */
					m3.build(&m2, s3, s3, s3);
					//skip building indices for m2 because they can be obtained directly
					//skip building indices for m3 because they can be obtained directly
					if (s3.size() > 4 * 28) {
						tbb::parallel_for(tbb::blocked_range<size_t>(0, s3.size(), 1), Loop2(ctx, s2, s3, m2, m3), tbb::auto_partitioner()); continue;
					}
					for (size_t i2_idx = 0; i2_idx < s3.size(); i2_idx++) { // loop-2 begin
						const IdType i2_id = s3[i2_idx];
						VertexSet m2_adj = m2.N(i2_idx);
						VertexSet m3_adj = m3.N(i2_idx);
						VertexSet s4 = m2_adj;
						if (s4.size() == 0) continue;
						/* VSet(4, 2) In-Edges: 0 1 2 Restricts: */
						VertexSet s5 = m3_adj.bounded(i2_id);
						if (s5.size() == 0) continue;
						/* VSet(5, 2) In-Edges: 0 1 2 Restricts: 0 1 2 */
						MiniGraphEager m4(false,false);
						/* Vertices = VSet(5) In-Edges: 0 1 2 Restricts: 0 1 2  | Intersect = VSet(4) In-Edges: 0 1 2 Restricts: */
						m4.build(&m2, s5, s4, s5);
						//skip building indices for m4 because they can be obtained directly
						if (s5.size() > 4 * 28) {
							tbb::parallel_for(tbb::blocked_range<size_t>(0, s5.size(), 1), Loop3(ctx, s4, s5, m4), tbb::auto_partitioner()); continue;
						}
						for (size_t i3_idx = 0; i3_idx < s5.size(); i3_idx++) { // loop-3 begin
							const IdType i3_id = s5[i3_idx];
							VertexSet m4_adj = m4.N(i3_idx);
							VertexSet s6 = m4_adj;
							if (s6.size() == 0) continue;
							/* VSet(6, 3) In-Edges: 0 1 2 3 Restricts: */
							VertexSet s7 = s5.intersect(m4_adj, m4_adj.vid());
							if (s7.size() == 0) continue;
							/* VSet(7, 3) In-Edges: 0 1 2 3 Restricts: 0 1 2 3 */
							auto m4_s7 = m4.indices(s7);
							for (size_t i4_idx = 0; i4_idx < s7.size(); i4_idx++) { // loop-4 begin
								VertexSet m4_adj = m4.N(m4_s7[i4_idx]);
								const IdType i4_id = s7[i4_idx];
								VertexSet s8 = s6.intersect(m4_adj);
								if (s8.size() == 0) continue;
								/* VSet(8, 4) In-Edges: 0 1 2 3 4 Restricts: */
								counter += 1ll * s8.size() * s8.size();
								/* Val: 1 | Group: (0),(1) | Comp: |VSet(8)|*|VSet(8)| */
								counter += -1ll * s8.size();
								/* Val: -1 | Group: (0 1) | Comp: |VSet(8) & VSet(8)| */
								} // loop-4 end
							} // loop-3 end
						} // loop-2 end
					} // loop-1 end
				} // loop-0 end
		} // operator end
	}; // Loop

	void plan(const GraphType* _graph, Context& ctx){ // plan 
		ctx.tick_begin = tbb::tick_count::now();
		ctx.iep_redundency = 2;
		graph = _graph;
		MiniGraphIF::DATA_GRAPH = graph;
		VertexSetType::MAX_DEGREE = graph->get_maxdeg();
		tbb::parallel_for(tbb::blocked_range<size_t>(0, graph->get_vnum()), Loop0(ctx), tbb::simple_partitioner());
	} // plan
} // minigraph
