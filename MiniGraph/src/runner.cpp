#include "codegen_output/plan.h"
#include "configure.h"
#include "common.h"
#include <filesystem>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>
#include <iostream>
#include <thread>
#include <chrono>
#include <mutex>
#include <math.h>
#include <condition_variable>
using namespace std::chrono_literals;
namespace minigraph {
    template<typename T>
    inline void read_file(std::filesystem::path path, T *&pointer, uint64_t num_elements) {
        CHECK(std::filesystem::is_regular_file(path)) << "File does not exists: " << path;
        if (pointer != nullptr) free(pointer);
        const size_t num_bytes = sizeof(T) * num_elements;
        pointer = new T[num_elements];
        std::ifstream file;
        file.open(path, std::ios::binary | std::ios::in);
        file.read(reinterpret_cast<char *>(pointer), num_bytes);
        CHECK(file.gcount() == (long int) num_bytes) << "Only read " << ToReadableSize(file.gcount()) << " out of "
                                                     << ToReadableSize(num_bytes) << " from " << path;
        file.close();
    };

    template<typename T>
    inline void mmap_file(std::filesystem::path path, T *&pointer, uint64_t num_elements) {
        CHECK(std::filesystem::is_regular_file(path)) << "File does not exists: " << path;
        if (pointer != nullptr) free(pointer);
        int fd = open(path.c_str(), O_RDONLY, 0);
        CHECK(fd != -1) << "Failed to open: " << path;
        pointer = static_cast<T *>(
                mmap(nullptr, sizeof(T) * num_elements, PROT_READ, MAP_SHARED, fd, 0));
        CHECK(pointer != MAP_FAILED) << "Failed to map file: " << path;
        CHECK(close(fd) == 0) << "Failed to close file: " << path;
    };

    inline GraphType *load_bin(std::string _in_dir, bool _mmap) {
        GraphType *out = new GraphType;
        MetaData m_meta;
        m_meta.read(_in_dir);

        out->m_mmap = _mmap;
        out->num_vertex = m_meta.num_vertex;
        out->num_edge = m_meta.num_edge;
        out->num_triangle = m_meta.num_triangle;
        out->max_offset = m_meta.max_offset;
        out->max_degree = m_meta.max_degree;
        out->max_triangle = m_meta.max_triangle;

        std::filesystem::path indicesFile = _in_dir;
        if (sizeof(IdType) == sizeof(uint64_t)) indicesFile /= Constant::kIndicesU64File;
        else if (sizeof(IdType) == sizeof(uint32_t)) indicesFile /= Constant::kIndicesU32File;
        else exit(-1 && "unsupported IdType");

        read_file<uint64_t>(std::filesystem::path{_in_dir} / Constant::kIndptrU64File, out->m_indptr,
                            m_meta.num_vertex + 1);
        read_file<uint64_t>(std::filesystem::path{_in_dir} / Constant::kOffsetU64File, out->m_offset,
                            m_meta.num_vertex);
        read_file<uint64_t>(std::filesystem::path{_in_dir} / Constant::kTriangleU64File, out->m_triangles,
                            m_meta.num_vertex);
        if (_mmap) {
            mmap_file<IdType>(indicesFile, out->m_indices, m_meta.num_edge);
        } else {
            read_file<IdType>(indicesFile, out->m_indices, m_meta.num_edge);
        }
        return out;
    }

    std::ostream &operator<<(std::ostream &os, const VertexSetType &dt) {
        if (dt.vid() == Constant::EmptyID<IdType>()) {
            os << "VertexSet(-1)\t=\t[";
        } else {
            os << "VertexSet(" << dt.vid() << ")\t=\t[";
        }
        constexpr uint64_t expected_to_show = 10;
        uint64_t actual_show = std::min(expected_to_show, dt.size());
        for (size_t i = 0; i < actual_show; i++) {
            os << dt[i] << ",";
        }
        if (expected_to_show < dt.size()) os << "...";
        os << "]\t";
        os << "(size=" << dt.size() << " pooled=" << dt.pooled() << ")\n";
        return os;
    }
}

void plan_2hrs(minigraph::GraphType* graph, minigraph::Context& ctx) {
    std::mutex m;
    std::condition_variable cv;
    std::thread t([&cv, &graph, &ctx](){
        minigraph::plan(graph, ctx);
        cv.notify_one();
    });
    t.detach();
    {
        std::unique_lock<std::mutex> l(m);
        // timeout after 1 day = 24 * 3600s = 86400s
        // timeout after 2 hours = 2 * 3600s = 7200s
        if (cv.wait_for(l, 7200s) == std::cv_status::timeout) throw std::runtime_error("Timeout");
    }
};

void plan_24hrs(minigraph::GraphType* graph, minigraph::Context& ctx) {
    std::mutex m;
    std::condition_variable cv;
    std::thread t([&cv, &graph, &ctx](){
        minigraph::plan(graph, ctx);
        cv.notify_one();
    });
    t.detach();
    {
        std::unique_lock<std::mutex> l(m);
        // timeout after 1 day = 24 * 3600s = 86400s
        // timeout after 2 hours = 2 * 3600s = 7200s
        if (cv.wait_for(l, 86400s) == std::cv_status::timeout) throw std::runtime_error("Timeout");
    }
};

int main(int argc, char *argv[]){
    using namespace minigraph;
    int expId = std::stoi(argv[1]);
    std::string in_dir{argv[2]};
    Timer t;
    const int processor_count = std::thread::hardware_concurrency();
    LOG(MSG) << "Threads=" << processor_count;
    GraphType *graph = load_bin(in_dir, false);
    LOG(MSG) << "LoadTime(s)=" << t.Passed();
    bool time_out = false;
    double seconds = 24 * 3600;
    Context ctx(processor_count);

    RunnerLog log;
    long long result{0};

    std::vector<int> skipIds = {93,138,139,140};
    for (int skipId: skipIds) {
        if (expId == skipId){
            log.finished = false;
            log.expId = expId;
            log.result = result;
            log.numThread = processor_count;
            log.vertexAllocated = VertexSetType::TOTAL_ALLOCATED;
            log.miniGraphAllocated = VertexSetType ::TOTAL_ALLOCATED;
            log.threadTimeSTD = 0.0;
            log.save(PROJECT_LOG_DIR);
            return 0;
        }
    }

    // Start Running Experiment (24-hours at most)
    t.Reset();
    try {
        plan_24hrs(graph, ctx);
    } catch (std::runtime_error& e) {
        time_out = true;
    }

    if (time_out) {
        result = ctx.get_result();
        seconds = t.Passed();

        log.finished = false;
        log.expId = expId;
        log.result = result;
        log.numThread = processor_count;
        log.vertexAllocated = VertexSetType::TOTAL_ALLOCATED;
        log.miniGraphAllocated = VertexSetType ::TOTAL_ALLOCATED;
        log.runTime = seconds;
        log.threadMinTime = seconds;
        log.threadMeanTime = seconds;
        log.threadMaxTime = seconds;
        log.threadTimeSTD = 0.0;

        LOG(MSG) << "CODE_EXECUTION_TIME(s)=" << "Timeout";
        LOG(MSG) << "RESULT=" << ctx.get_result();
        LOG(MSG) << "Throughput=" << result / seconds;
        LOG(MSG) << "VertexSetAllocated=" << ToReadableSize(VertexSetType::TOTAL_ALLOCATED);
        LOG(MSG) << "MiniGraphAllocated=" << ToReadableSize(MiniGraphPool::TOTAL_ALLOCATED);
    } else {
        result = ctx.get_result();
        seconds = t.Passed();

        log.finished = true;
        log.expId = expId;
        log.result = result;
        log.numThread = processor_count;
        log.vertexAllocated = VertexSetType::TOTAL_ALLOCATED;
        log.miniGraphAllocated = MiniGraphPool::TOTAL_ALLOCATED;
        log.runTime = seconds;
        log.threadMinTime = ctx.get_min_time();
        log.threadMeanTime = ctx.get_mean_time();
        log.threadMaxTime = ctx.get_max_time();
        log.threadTimeSTD = sqrt(ctx.get_var_time());

        LOG(MSG) << "CODE_EXECUTION_TIME(s)=" << seconds;
        LOG(MSG) << "RESULT=" << ctx.get_result();
        LOG(MSG) << "Throughput=" << ctx.get_result() / seconds;
        // LOG(MSG) << "ThreadMeanTime=" << ctx.get_mean_time() << "s";
        // LOG(MSG) << "ThreadMinTime=" << ctx.get_min_time() << "s";
        // LOG(MSG) << "ThreadMaxTime=" << ctx.get_max_time() << "s";
        // LOG(MSG) << "TimeSTD=" << sqrt(ctx.get_var_time());
        LOG(MSG) << "VertexSetAllocated=" << ToReadableSize(VertexSetType::TOTAL_ALLOCATED);
        LOG(MSG) << "MiniGraphAllocated=" << ToReadableSize(MiniGraphPool::TOTAL_ALLOCATED);
    }
    log.save(PROJECT_LOG_DIR);
}