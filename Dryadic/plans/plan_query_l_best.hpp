const int n_counters = 1;
std::vector<std::vector<uint64_t>> global_counters(omp_get_max_threads());
const vidType vertices = g.V();
#pragma omp parallel
{
  std::vector<uint64_t> &counter = global_counters.at(omp_get_thread_num());
  counter.resize(n_counters, 0);
    VertexSet l0 = g.L(0);
    #pragma omp for schedule(dynamic,64) nowait
    for(vidType v0:l0)    {
      VertexSet y0l0 = g.N(v0,0);
      VertexSet y0f0l0 = bounded(y0l0, v0);
      for(vidType v1:y0f0l0)      {
        VertexSet y1l0 = g.N(v1,0);
        VertexSet y0y1l0 = intersection_set(y0l0, y1l0);
        VertexSet y0f0y1f1l0 = intersection_set(y0f0l0, y1l0, v1);
        for(vidType v2:y0f0y1f1l0)        {
          VertexSet y2l0 = g.N(v2,0);
          VertexSet y0y1y2l0 = intersection_set(y0y1l0, y2l0);
          VertexSet y0f0y1f1y2f2l0 = intersection_set(y0f0y1f1l0, y2l0, v2);
          for(vidType v3:y0f0y1f1y2f2l0)          {
            VertexSet y3l0 = g.N(v3,0);
            VertexSet y0y1y2y3l0 = intersection_set(y0y1y2l0, y3l0);
            VertexSet y0f0y1f1y2f2y3f3l0 = intersection_set(y0f0y1f1y2f2l0, y3l0, v3);
            for(vidType v4:y0f0y1f1y2f2y3f3l0)            {
              VertexSet y4l0 = g.N(v4,0);
              VertexSet y0y1y2y3y4l0 = intersection_set(y0y1y2y3l0, y4l0);
              for(vidType v5:y0y1y2y3y4l0)              {
                VertexSet y5l0 = g.N(v5,0);
                counter[0] += difference_num(y0y1y2y3y4l0, y5l0, v5);
              }
            }
          }
        }
      }
    }
}//217
const uint64_t data_complexity = 1.44948e+07;
const uint64_t time_complexity = 1.93507e+07;
