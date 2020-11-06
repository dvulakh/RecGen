
/********************************************************************
* Defines some general macros for writing to log files
********************************************************************/

#ifndef LOGGING_H
#define LOGGING_H

#include <cstdio>
#include <string>
#include <chrono>

// Accessing the settings bits of a logging object
#define IS(S) (settings & (S))
#define LOG_WORK 1LL << 0
#define LOG_DATA 1LL << 1
#define VER_WORK 1LL << 2
#define VER_DATA 1LL << 3

// Logging macros
#define TPLUS(t0) (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - (t0)) / 1000000000.0)
#define PRINT_HEADER(s) (std::string(30, '=') + " " + s + " " + std::string(30, '='))
#define GPRINTF(s, ...) (std::string("[%f]\t") + s + std::string("\n")).c_str(), TPLUS(this->start_time), __VA_ARGS__
#define GPRINT(s) (std::string("[%f]\t") + s + std::string("\n")).c_str(), TPLUS(this->start_time)
#define FPRINTF(f, s, v, ...) { if(v) { std::fprintf(f, GPRINTF(s, __VA_ARGS__)); std::fflush(f); } }
#define FPRINT(f, s, v) { if(v) { std::fprintf(f, GPRINT(s)); std::fflush(f); } }
#define PRINTF(s, v, ...) { if(v) { std::printf(GPRINTF(s, __VA_ARGS__)); std::fflush(stdout); } }
#define PRINT(s, v) { if(v) { std::printf(GPRINT(s)); std::fflush(stdout); } }
#define WPRINTF(s, ...) { FPRINTF(this->work_log, s, IS(LOG_WORK), __VA_ARGS__) PRINTF(s, IS(VER_WORK), __VA_ARGS__) }
#define DPRINTF(s, ...) { FPRINTF(this->data_log, s, IS(LOG_DATA), __VA_ARGS__) PRINTF(s, IS(VER_DATA), __VA_ARGS__) }
#define WPRINT(s) { FPRINT(this->work_log, s, IS(LOG_WORK)) PRINT(s, IS(VER_WORK)) }
#define DPRINT(s) { FPRINT(this->data_log, s, IS(LOG_DATA)) PRINT(s, IS(VER_DATA)) }

// Include the necessary members
#define MAKE_LOGGABLE \
protected: \
std::string work_path, data_path; \
std::FILE *work_log, *data_log; \
std::chrono::high_resolution_clock::time_point start_time; \
public: \
void set_work_path(std::string work_path) { this->work_path = work_path; } \
void set_data_path(std::string data_path) { this->data_path = data_path; } \
long long settings;

// Logging command-line flags
#define LOG_FLAG_READ(flg_rd, logger) \
fr.add_flag("verbose", 'v', 0, [&](std::vector<std::string> v, void* p) { logger->settings |= VER_WORK | VER_DATA; }); \
fr.add_flag("dataf", 'D', 1, [&](std::vector<std::string> v, void* p) { logger->set_data_path(v[0]); }); \
fr.add_flag("workf", 'W', 1, [&](std::vector<std::string> v, void* p) { logger->set_work_path(v[0]); }); \
fr.add_flag("data", 'd', 0, [&](std::vector<std::string> v, void* p) { logger->settings |= LOG_DATA | VER_DATA, logger->settings &= ~(LOG_WORK | VER_WORK); }); \
fr.add_flag("work", 'w', 0, [&](std::vector<std::string> v, void* p) { logger->settings |= LOG_WORK | VER_WORK, logger->settings &= ~(LOG_DATA | VER_DATA); }); \
fr.add_flag("silent", 's', 0, [&](std::vector<std::string> v, void* p) { logger->settings &= ~(LOG_DATA | LOG_WORK | VER_DATA | VER_WORK); }); \

#endif
