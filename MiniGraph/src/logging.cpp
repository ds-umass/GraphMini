//
// Created by Ubuntu on 1/10/23.
//

#include "logging.h"
#include "constant.h"
#include <cstdio>
#include <string>
#include <iostream>
#include <chrono>
#include <iomanip>
namespace minigraph {
    std::string ToReadableSize(size_t nbytes) {
        char buf[Constant::kBufferSize];
        if (nbytes >= Constant::kGigabytes) {
            double new_size = (float)nbytes / Constant::kGigabytes;
            snprintf(buf, Constant::kBufferSize, "%.2lf GB", new_size);
            return std::string(buf);
        } else if (nbytes >= Constant::kMegabytes) {
            double new_size = (float)nbytes / Constant::kMegabytes;
            snprintf(buf, Constant::kBufferSize, "%.2lf MB", new_size);
            return std::string(buf);
        } else if (nbytes >= Constant::kKilobytes) {
            double new_size = (float)nbytes / Constant::kKilobytes;
            snprintf(buf, Constant::kBufferSize, "%.2lf KB", new_size);
            return std::string(buf);
        } else {
            double new_size = (float)nbytes;
            snprintf(buf, Constant::kBufferSize, "%.2lf Bytes", new_size);
            return std::string(buf);
        }
    };
}

namespace minigraph {
    std::string LOG_LEVELS[] = {"trace", "debug", "info", "msg", "warning", "error", "fatal"};

    LogMessage::LogMessage(const char* fname, int line, LogLevel severity)
            : fname_(fname), line_(line), severity_(severity) {}

    void LogMessage::GenerateLogMessage(bool log_time) {
        bool use_cout =
                static_cast<int>(severity_) <= static_cast<int>(LogLevel::MSG);
        std::ostream& os = use_cout ? std::cout : std::cerr;
        if (severity_ == LogLevel::MSG) {
            os << str() << std::endl;
        } else if (log_time) {
            auto now = std::chrono::system_clock::now();
            auto as_time_t = std::chrono::system_clock::to_time_t(now);

            auto duration = now.time_since_epoch();
            auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
            auto micros_remainder =
                    std::chrono::duration_cast<std::chrono::microseconds>(duration -
                                                                          seconds);

            const size_t time_buffer_size = 30;
            char time_buffer[time_buffer_size];
            strftime(time_buffer, time_buffer_size, "%Y-%m-%d %H:%M:%S",
                     localtime(&as_time_t));
            os << "[" << time_buffer << "." << std::setw(3) << micros_remainder.count()
               << ": "<< fname_ << ":" << line_ << "] " << "(" << LOG_LEVELS[static_cast<int>(severity_)] << ") " << str() << std::endl;
        } else {
            os << "[" << " " << fname_ << ":"
               << line_ << "] " "(" << LOG_LEVELS[static_cast<int>(severity_)] << ") " << str() << std::endl;
        }
    }

    LogMessage::~LogMessage() {
        static LogLevel min_log_level = MinLogLevelFromEnv();
        static bool log_time = LogTimeFromEnv();
        if (severity_ >= min_log_level) {
            GenerateLogMessage(log_time);
        }
    }

    LogMessageFatal::LogMessageFatal(const char* file, int line)
            : LogMessage(file, line, LogLevel::FATAL) {}

    LogMessageFatal::~LogMessageFatal() {
        static bool log_time = LogTimeFromEnv();
        GenerateLogMessage(log_time);
        exit(-1);
    }

    LogLevel ParseLogLevelStr(const char* env_var_val) {
        std::string min_log_level(env_var_val);
        if (min_log_level == "trace") {
            return LogLevel::TRACE;
        } else if (min_log_level == "debug") {
            return LogLevel::DEBUG;
        } else if (min_log_level == "info") {
            return LogLevel::INFO;
        } else if (min_log_level == "warning") {
            return LogLevel::WARNING;
        } else if (min_log_level == "error") {
            return LogLevel::ERROR;
        } else if (min_log_level == "fatal") {
            return LogLevel::FATAL;
        } else {
            return LogLevel::WARNING;
        }
    }

    LogLevel MinLogLevelFromEnv() {
        const char* env_var_val = getenv("MINIGRAPH_LOG_LEVEL");
        if (env_var_val == nullptr) {
            return LogLevel::INFO;
        }
        return ParseLogLevelStr(env_var_val);
    }

    bool LogTimeFromEnv() {
        const char* env_var_val = getenv("MINIGRAPH_LOG_HIDE_TIME");
        if (env_var_val != nullptr && std::strtol(env_var_val, nullptr, 10) > 0) {
            return false;
        } else {
            return true;
        }
    }
}