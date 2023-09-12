#include "plan_profile.h"
#include "oneapi/tbb/parallel_for.h"
namespace minigraph {
	uint64_t pattern_size() {return 6;}
	using namespace oneapi::tbb;
	static const Graph * graph;
	class Loop3
	{
	private:
		Context& ctx;
		// Parent Intermediates
		VertexSet& s5;
		// Iterate Set
		VertexSet& s6;
	public:
		Loop3(Context& _ctx, VertexSet& _s5, VertexSet& _s6):ctx{_ctx}, s5{_s5}, s6{_s6} {};
		void operator()(const blocked_range<size_t> &r) const {// operator begin
			const int worker_id = tbb::this_task_arena::current_thread_index();
			cc& counter = ctx.per_thread_result.at(worker_id);
			for (size_t i3_idx = r.begin(); i3_idx < r.end(); i3_idx++) { // loop-3begin
				ctx.profiler->set_cur_loop(3);
							const IdType i3_id = s6[i3_idx];
							VertexSet i3_adj = graph->N(i3_id);
				VertexSet s7 = s5.subtract(i3_adj);
				/* VSet(7, 3) In-Edges: 1 2 Restricts: */
				VertexSet s8 = s6.intersect(i3_adj, i3_adj.vid());
				/* VSet(8, 3) In-Edges: 0 1 3 Restricts: 3 */
				for (size_t i4_idx = 0; i4_idx < s8.size(); i4_idx++) { // loop-4 begin
					ctx.profiler->set_cur_loop(4);
								const IdType i4_id = s8[i4_idx];
								VertexSet i4_adj = graph->N(i4_id);
					counter += s7.subtract_cnt(i4_adj);
					/* VSet(9, 4) In-Edges: 1 2 Restricts: */
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
		VertexSet& s4;
		// Iterate Set
		VertexSet& s3;
	public:
		Loop2(Context& _ctx, VertexSet& _s2, VertexSet& _s4, VertexSet& _s3):ctx{_ctx}, s2{_s2}, s4{_s4}, s3{_s3} {};
		void operator()(const blocked_range<size_t> &r) const {// operator begin
			const int worker_id = tbb::this_task_arena::current_thread_index();
			cc& counter = ctx.per_thread_result.at(worker_id);
			for (size_t i2_idx = r.begin(); i2_idx < r.end(); i2_idx++) { // loop-2begin
				ctx.profiler->set_cur_loop(2);
						const IdType i2_id = s3[i2_idx];
						VertexSet i2_adj = graph->N(i2_id);
				VertexSet s5 = s2.intersect(i2_adj);
				/* VSet(5, 2) In-Edges: 1 2 Restricts: */
				VertexSet s6 = s4.subtract(i2_adj);
				/* VSet(6, 2) In-Edges: 0 1 Restricts: */
				if (s6.size() > 4 * 17) {
					parallel_for(blocked_range<size_t>(0, s6.size(), 1), Loop3(ctx, s5, s6), auto_partitioner()); continue;
				}
				for (size_t i3_idx = 0; i3_idx < s6.size(); i3_idx++) { // loop-3 begin
					ctx.profiler->set_cur_loop(3);
							const IdType i3_id = s6[i3_idx];
							VertexSet i3_adj = graph->N(i3_id);
					VertexSet s7 = s5.subtract(i3_adj);
					/* VSet(7, 3) In-Edges: 1 2 Restricts: */
					VertexSet s8 = s6.intersect(i3_adj, i3_adj.vid());
					/* VSet(8, 3) In-Edges: 0 1 3 Restricts: 3 */
					for (size_t i4_idx = 0; i4_idx < s8.size(); i4_idx++) { // loop-4 begin
						ctx.profiler->set_cur_loop(4);
								const IdType i4_id = s8[i4_idx];
								VertexSet i4_adj = graph->N(i4_id);
						counter += s7.subtract_cnt(i4_adj);
						/* VSet(9, 4) In-Edges: 1 2 Restricts: */
						} // loop-4 end
					} // loop-3 end
				} // loop-2 end
		} // operator end
	}; // Loop

	class Loop1
	{
	private:
		Context& ctx;
		// Adjacent Lists
		VertexSet& i0_adj;
		// Parent Intermediates
		VertexSet& s0;
		// Iterate Set
		VertexSet& s1;
	public:
		Loop1(Context& _ctx, VertexSet& _i0_adj, VertexSet& _s0, VertexSet& _s1):ctx{_ctx}, i0_adj{_i0_adj}, s0{_s0}, s1{_s1} {};
		void operator()(const blocked_range<size_t> &r) const {// operator begin
			const int worker_id = tbb::this_task_arena::current_thread_index();
			cc& counter = ctx.per_thread_result.at(worker_id);
			for (size_t i1_idx = r.begin(); i1_idx < r.end(); i1_idx++) { // loop-1begin
				ctx.profiler->set_cur_loop(1);
					const IdType i1_id = s1[i1_idx];
					VertexSet i1_adj = graph->N(i1_id);
				VertexSet s2 = i1_adj.subtract(i0_adj);
					if (s2.size() == 0) continue;
				/* VSet(2, 1) In-Edges: 1 Restricts: */
				VertexSet s3 = s0.subtract(i1_adj);
				/* VSet(3, 1) In-Edges: 0 Restricts: */
				VertexSet s4 = s0.intersect(i1_adj);
				/* VSet(4, 1) In-Edges: 0 1 Restricts: */
				if (s3.size() > 4 * 17) {
					parallel_for(blocked_range<size_t>(0, s3.size(), 1), Loop2(ctx, s2, s4, s3), auto_partitioner()); continue;
				}
				for (size_t i2_idx = 0; i2_idx < s3.size(); i2_idx++) { // loop-2 begin
					ctx.profiler->set_cur_loop(2);
						const IdType i2_id = s3[i2_idx];
						VertexSet i2_adj = graph->N(i2_id);
					VertexSet s5 = s2.intersect(i2_adj);
					/* VSet(5, 2) In-Edges: 1 2 Restricts: */
					VertexSet s6 = s4.subtract(i2_adj);
					/* VSet(6, 2) In-Edges: 0 1 Restricts: */
					if (s6.size() > 4 * 17) {
						parallel_for(blocked_range<size_t>(0, s6.size(), 1), Loop3(ctx, s5, s6), auto_partitioner()); continue;
					}
					for (size_t i3_idx = 0; i3_idx < s6.size(); i3_idx++) { // loop-3 begin
						ctx.profiler->set_cur_loop(3);
							const IdType i3_id = s6[i3_idx];
							VertexSet i3_adj = graph->N(i3_id);
						VertexSet s7 = s5.subtract(i3_adj);
						/* VSet(7, 3) In-Edges: 1 2 Restricts: */
						VertexSet s8 = s6.intersect(i3_adj, i3_adj.vid());
						/* VSet(8, 3) In-Edges: 0 1 3 Restricts: 3 */
						for (size_t i4_idx = 0; i4_idx < s8.size(); i4_idx++) { // loop-4 begin
							ctx.profiler->set_cur_loop(4);
								const IdType i4_id = s8[i4_idx];
								VertexSet i4_adj = graph->N(i4_id);
							counter += s7.subtract_cnt(i4_adj);
							/* VSet(9, 4) In-Edges: 1 2 Restricts: */
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
		void operator()(const blocked_range<size_t> &r) const {// operator begin
			const int worker_id = tbb::this_task_arena::current_thread_index();
			cc& counter = ctx.per_thread_result.at(worker_id);
			for (size_t i0_id = r.begin(); i0_id < r.end(); i0_id++) { // loop-0begin
				ctx.profiler->set_cur_loop(0);
				VertexSet i0_adj = graph->N(i0_id);
				VertexSet s0 = i0_adj;
				if (s0.size() == 0) continue;
				/* VSet(0, 0) In-Edges: 0 Restricts: */
				VertexSet s1 = s0.bounded(i0_id);
				/* VSet(1, 0) In-Edges: 0 Restricts: 0 */
				if (s1.size() > 4 * 17) {
					parallel_for(blocked_range<size_t>(0, s1.size(), 1), Loop1(ctx, i0_adj, s0, s1), auto_partitioner()); continue;
				}
				for (size_t i1_idx = 0; i1_idx < s1.size(); i1_idx++) { // loop-1 begin
					ctx.profiler->set_cur_loop(1);
					const IdType i1_id = s1[i1_idx];
					VertexSet i1_adj = graph->N(i1_id);
					VertexSet s2 = i1_adj.subtract(i0_adj);
					if (s2.size() == 0) continue;
					/* VSet(2, 1) In-Edges: 1 Restricts: */
					VertexSet s3 = s0.subtract(i1_adj);
					/* VSet(3, 1) In-Edges: 0 Restricts: */
					VertexSet s4 = s0.intersect(i1_adj);
					/* VSet(4, 1) In-Edges: 0 1 Restricts: */
					if (s3.size() > 4 * 17) {
						parallel_for(blocked_range<size_t>(0, s3.size(), 1), Loop2(ctx, s2, s4, s3), auto_partitioner()); continue;
					}
					for (size_t i2_idx = 0; i2_idx < s3.size(); i2_idx++) { // loop-2 begin
						ctx.profiler->set_cur_loop(2);
						const IdType i2_id = s3[i2_idx];
						VertexSet i2_adj = graph->N(i2_id);
						VertexSet s5 = s2.intersect(i2_adj);
						/* VSet(5, 2) In-Edges: 1 2 Restricts: */
						VertexSet s6 = s4.subtract(i2_adj);
						/* VSet(6, 2) In-Edges: 0 1 Restricts: */
						if (s6.size() > 4 * 17) {
							parallel_for(blocked_range<size_t>(0, s6.size(), 1), Loop3(ctx, s5, s6), auto_partitioner()); continue;
						}
						for (size_t i3_idx = 0; i3_idx < s6.size(); i3_idx++) { // loop-3 begin
							ctx.profiler->set_cur_loop(3);
							const IdType i3_id = s6[i3_idx];
							VertexSet i3_adj = graph->N(i3_id);
							VertexSet s7 = s5.subtract(i3_adj);
							/* VSet(7, 3) In-Edges: 1 2 Restricts: */
							VertexSet s8 = s6.intersect(i3_adj, i3_adj.vid());
							/* VSet(8, 3) In-Edges: 0 1 3 Restricts: 3 */
							for (size_t i4_idx = 0; i4_idx < s8.size(); i4_idx++) { // loop-4 begin
								ctx.profiler->set_cur_loop(4);
								const IdType i4_id = s8[i4_idx];
								VertexSet i4_adj = graph->N(i4_id);
								counter += s7.subtract_cnt(i4_adj);
								/* VSet(9, 4) In-Edges: 1 2 Restricts: */
								} // loop-4 end
							} // loop-3 end
						} // loop-2 end
					} // loop-1 end
				} // loop-0 end
		} // operator end
	}; // Loop

	void plan(const GraphType* _graph, Context& ctx){ // plan 
		VertexSet::profiler = ctx.profiler;
		ctx.tick_begin = tick_count::now();
		ctx.iep_redundency = 0;
		graph = _graph;
		VertexSetType::MAX_DEGREE = graph->get_maxdeg();
		parallel_for(blocked_range<size_t>(0, graph->get_vnum()), Loop0(ctx), simple_partitioner());
	} // plan
} // minigraph
