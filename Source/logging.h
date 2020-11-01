
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
std::string work_path, data_path; \
std::FILE *work_log, *data_log; \
std::chrono::high_resolution_clock::time_point start_time; \
long long settings;

#endif
