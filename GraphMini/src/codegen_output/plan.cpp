#include "plan.h"
namespace minigraph {
uint64_t pattern_size() { return 7; }
void plan(const GraphType *graph, Context &ctx) {
  VertexSetType::MAX_DEGREE = graph->get_maxdeg();
#pragma omp parallel num_threads(ctx.num_threads) default(none)                \
    shared(ctx, graph)
  { // pragma parallel
    cc &counter = ctx.per_thread_result.at(omp_get_thread_num());
    cc &handled = ctx.per_thread_handled.at(omp_get_thread_num());
    double start = omp_get_wtime();
    ctx.iep_redundency = 0;
#pragma omp for schedule(dynamic, 1) nowait
    for (IdType i0_id = 0; i0_id < graph->get_vnum(); i0_id++) { // loop-0 begin
      VertexSet i0_adj = graph->N(i0_id);
      VertexSet s0 = i0_adj;
      if (s0.size() == 0)
        continue;
      /* VSet(0, 0) In-Edges: 0 Restricts: */
      VertexSet s1 = s0.bounded(i0_id);
      /* VSet(1, 0) In-Edges: 0 Restricts: 0 */
      for (size_t i1_idx = 0; i1_idx < s1.size(); i1_idx++) { // loop-1 begin
        const IdType i1_id = s1[i1_idx];
        VertexSet i1_adj = graph->N(i1_id);
        VertexSet s2 = s0.intersect(i1_adj);
        /* VSet(2, 1) In-Edges: 0 1 Restricts: */
        VertexSet s3 = s1.intersect(i1_adj, i1_adj.vid());
        /* VSet(3, 1) In-Edges: 0 1 Restricts: 0 1 */
        for (size_t i2_idx = 0; i2_idx < s3.size(); i2_idx++) { // loop-2 begin
          const IdType i2_id = s3[i2_idx];
          VertexSet i2_adj = graph->N(i2_id);
          VertexSet s4 = s2.intersect(i2_adj);
          /* VSet(4, 2) In-Edges: 0 1 2 Restricts: */
          VertexSet s5 = s3.intersect(i2_adj, i2_adj.vid());
          /* VSet(5, 2) In-Edges: 0 1 2 Restricts: 0 1 2 */
          for (size_t i3_idx = 0; i3_idx < s5.size();
               i3_idx++) { // loop-3 begin
            const IdType i3_id = s5[i3_idx];
            VertexSet i3_adj = graph->N(i3_id);
            VertexSet s6 = s4.intersect(i3_adj);
            /* VSet(6, 3) In-Edges: 0 1 2 3 Restricts: */
            VertexSet s7 = s5.intersect(i3_adj, i3_adj.vid());
            /* VSet(7, 3) In-Edges: 0 1 2 3 Restricts: 0 1 2 3 */
            for (size_t i4_idx = 0; i4_idx < s7.size();
                 i4_idx++) { // loop-4 begin
              const IdType i4_id = s7[i4_idx];
              VertexSet i4_adj = graph->N(i4_id);
              VertexSet s8 = s6.intersect(i4_adj);
              /* VSet(8, 4) In-Edges: 0 1 2 3 4 Restricts: */
              for (size_t i5_idx = 0; i5_idx < s8.size();
                   i5_idx++) { // loop-5 begin
                const IdType i5_id = s8[i5_idx];
                VertexSet i5_adj = graph->N(i5_id);
                counter += s8.subtract_cnt(i5_adj, i5_adj.vid());
                /* VSet(9, 5) In-Edges: 0 1 2 3 4 Restricts: 5 */
              } // loop-5 end
            }   // loop-4 end
          }     // loop-3 end
        }       // loop-2 end
      }         // loop-1 end
      handled += 1;
    } // loop-0 end
    ctx.per_thread_time.at(omp_get_thread_num()) = omp_get_wtime() - start;
  } // pragma parallel
} // plan
} // namespace minigraph
extern "C" void plan(const minigraph::GraphType *graph,
                     minigraph::Context &ctx) {
  return minigraph::plan(graph, ctx);
};