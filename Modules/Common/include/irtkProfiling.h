/* The Image Registration Toolkit (IRTK)
 *
 * Copyright 2008-2015 Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. */

#ifndef _IRTKPROFILING_H
#define _IRTKPROFILING_H

#include <iostream>
#include <sstream>
#include <time.h>

#ifdef HAS_TBB
#  include <tbb/tick_count.h>
#endif


// -----------------------------------------------------------------------------
enum TimeUnit {TIME_IN_DEFAULT_UNIT, TIME_IN_MILLISECONDS, TIME_IN_SECONDS};

// =============================================================================
// Global profiling options
// =============================================================================

/// Enable/disable profiling of execution time.
///
/// Should be set in the main function before any processing starts, e.g.,
/// depending on a command-line flag (for example -v -verbose).
/// If less or equal to zero, no timing measurements are printed to screen.
/// Otherwise, whether a timing measure is output or not depends on the
/// set debugging level.
extern int debug_time;

/// Time unit to use for output of time measurements.
extern TimeUnit debug_time_unit;

// =============================================================================
// Command help
// =============================================================================

/// Check if given option is a profiling option
bool IsProfilingOption(const char *);

/// Parse profiling option
void ParseProfilingOption(int &, int &, char *[]);

/// Print profiling command-line options
void PrintProfilingOptions(std::ostream &);

// =============================================================================
// CPU Profiling
// =============================================================================

// -----------------------------------------------------------------------------
/// Print elapsed time for given section
void PrintElapsedTime(const char *, double, TimeUnit = TIME_IN_SECONDS);

// -----------------------------------------------------------------------------
/// Start measurement of execution time of current code block
///
/// @code
/// {
///   IRTK_START_TIMING();
///   // do some work here
///   IRTK_END_TIMING("example section");
/// }
/// @endcode
///
/// @sa IRTK_END_TIMING
#ifdef USE_TIMING
#  ifdef HAS_TBB
#    define IRTK_START_TIMING()   tbb::tick_count t_start = tbb::tick_count::now()
#  else
#    define IRTK_START_TIMING()   clock_t t_start = clock()
#  endif
#else
#  define IRTK_START_TIMING()   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// Reset measurement of starting execution time of current code block
///
/// @code
/// {
///   IRTK_START_TIMING();
///   // do some work here
///   IRTK_END_TIMING("first part");
///   IRTK_RESET_TIMING();
///   // do some work here
///   IRTK_END_TIMING("second part");
/// }
/// @endcode
///
/// @sa IRTK_END_TIMING
#ifdef USE_TIMING
#  ifdef HAS_TBB
#    define IRTK_RESET_TIMING()   t_start = tbb::tick_count::now()
#  else
#    define IRTK_RESET_TIMING()   t_start = clock()
#  endif
#else
#  define IRTK_RESET_TIMING()   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// End measurement of execution time of current code block.
///
/// @code
/// {
///   IRTK_START_TIMING();
///   // do some work here
///   IRTK_END_TIMING("example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at compile time depending on the
///       USE_TIMING flag.
///
/// @sa IRTK_START_TIMING
#ifdef USE_TIMING
#  ifdef HAS_TBB
#    define IRTK_END_TIMING(section)                                           \
       do {                                                                    \
         std::ostringstream oss;                                               \
         oss << section;                                                       \
         PrintElapsedTime(oss.str().c_str(),                                   \
                          (tbb::tick_count::now() - t_start).seconds());       \
       } while (false)
#  else
#    define IRTK_END_TIMING(section)                                           \
       do {                                                                    \
         std::ostringstream oss;                                               \
         oss << section;                                                       \
         PrintElapsedTime(oss.str().c_str(),                                   \
                          static_cast<double>(clock() - t_start)               \
                        / static_cast<double>(CLOCKS_PER_SEC));                \
       } while (false)
#  endif
#else
#  define IRTK_END_TIMING(section)   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// End measurement of execution time of current code block.
///
/// In the following example, the timing measurement are only printed if
/// the global execution time debugging level is greater or equal the
/// specified debugging level. Hence, if debug_time is 1, only the time of
/// the following example section is printed. If debug_time is greater than 1,
/// also the times needed for each execution of inner block are printed together
/// with the summarizing total execution time of the example section. Otherwise,
/// if debug_time is less than 1, no time measurements are printed at all.
/// @code
/// {
///   IRTK_START_TIMING();
///   // do some work here
///   {
///     IRTK_START_TIMING();
///     // do some part of the work here
///     IRTK_DEBUG_TIMING(2, "example part");
///   }
///   // possibly combine results and maybe do some clean up work here
///   IRTK_DEBUG_TIMING(1, "example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at runtime depending on the global
///       variable debug_time.
///
/// @sa IRTK_START_TIMING
#ifdef USE_TIMING
#  ifdef HAS_TBB
#    define IRTK_DEBUG_TIMING(level, section)                                  \
       do {                                                                    \
         if (debug_time >= level) {                                            \
           std::ostringstream oss;                                             \
           oss << section;                                                     \
           PrintElapsedTime(oss.str().c_str(),                                 \
                            (tbb::tick_count::now() - t_start).seconds());     \
         }                                                                     \
       } while (false)
#  else
#    define IRTK_DEBUG_TIMING(level, section)                                  \
       do {                                                                    \
         if (debug_time >= level) {                                            \
           std::ostringstream oss;                                             \
           oss << section;                                                     \
           PrintElapsedTime(oss.str().c_str(),                                 \
                            static_cast<double>(clock() - t_start)             \
                          / static_cast<double>(CLOCKS_PER_SEC));              \
         }                                                                     \
       } while (false)
#  endif
#else
#  define IRTK_DEBUG_TIMING(level, section)   do {} while (false)
#endif

// =============================================================================
// GPU Profiling
// =============================================================================

// -----------------------------------------------------------------------------
/// Start measurement of execution time of current code block.
///
/// @code
/// {
///   IRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   IRTKCU_END_TIMING("example section");
/// }
/// @endcode
///
/// @sa IRTKCU_END_TIMING
#ifdef USE_TIMING
#  define IRTKCU_START_TIMING()                                                \
           cudaEvent_t e_start, e_stop;                                        \
           CudaSafeCall( cudaEventCreate(&e_start) );                          \
           CudaSafeCall( cudaEventCreate(&e_stop) );                           \
           CudaSafeCall( cudaEventRecord(e_start, 0) )
#else
#  define IRTKCU_START_TIMING()   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// Reset measurement of starting execution time of current code block
///
/// @code
/// {
///   IRTKCU_START_TIMING();
///   // do some work here
///   IRTKCU_END_TIMING("first part");
///   IRTKCU_RESET_TIMING();
///   // do some work here
///   IRTKCU_END_TIMING("second part");
/// }
/// @endcode
///
/// @sa IRTKCU_END_TIMING
#ifdef USE_TIMING
#  define IRTKCU_RESET_TIMING()   CudaSafeCall( cudaEventRecord(e_start, 0) )
#else
#  define IRTKCU_RESET_TIMING()   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// Print interim measurement of execution time of current code block.
///
/// @code
/// {
///   IRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   IRTKCU_INTERIM_TIMING("example kernel");
///   CudaSafeCall( cudaDeviceSynchronize() );
///   IRTKCU_END_TIMING("example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at compile time depending on the
///       USE_TIMING flag.
///
/// @sa IRTKCU_START_TIMING
#ifdef USE_TIMING
#  define IRTKCU_INTERIM_TIMING(section)                                       \
     do {                                                                      \
       float t_elapsed;                                                        \
       CudaSafeCall( cudaEventRecord(e_stop, 0) );                             \
       CudaSafeCall( cudaEventSynchronize(e_stop) );                           \
       CudaSafeCall( cudaEventElapsedTime(&t_elapsed, e_start, e_stop) );      \
       std::ostringstream oss;                                                 \
       oss << section << " [interim]";                                         \
       PrintElapsedTime(oss.str().c_str(), t_elapsed, TIME_IN_MILLISECONDS);   \
     } while (false)
#else
#  define IRTKCU_INTERIM_TIMING(section)   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// End measurement of execution time of current code block.
///
/// @code
/// {
///   IRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   IRTKCU_END_TIMING("example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at compile time depending on the
///       USE_TIMING flag.
///
/// @sa IRTKCU_START_TIMING
#ifdef USE_TIMING
#  define IRTKCU_END_TIMING(section)                                           \
     do {                                                                      \
       float t_elapsed;                                                        \
       CudaSafeCall( cudaEventRecord(e_stop, 0) );                             \
       CudaSafeCall( cudaEventSynchronize(e_stop) );                           \
       CudaSafeCall( cudaEventElapsedTime(&t_elapsed, e_start, e_stop) );      \
       CudaSafeCall( cudaEventDestroy(e_start) );                              \
       CudaSafeCall( cudaEventDestroy(e_stop) );                               \
       std::ostringstream oss;                                                 \
       oss << section << " [GPU]";                                             \
       PrintElapsedTime(oss.str().c_str(), t_elapsed, TIME_IN_MILLISECONDS);   \
     } while (false)
#else
#  define IRTKCU_END_TIMING(section)   do {} while (false)
#endif

// -----------------------------------------------------------------------------
/// Print interim measurement of execution time of current code block.
///
/// In the following example, the timing measurement are only printed if
/// the global execution time debugging level is greater or equal the
/// specified debugging level. Hence, if debug_time is 1, only the time of
/// the following example section is printed. If debug_time is greater than 1,
/// also the times needed for each execution of inner block are printed together
/// with the summarizing total execution time of the example section. Otherwise,
/// if debug_time is less than 1, no time measurements are printed at all.
/// @code
/// {
///   IRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   IRTKCU_DEBUG_INTERIM_TIMING(1, "kernel name");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at runtime depending on the global
///       variable debug_time.
///
/// @sa IRTKCU_START_TIMING
#ifdef USE_TIMING
#  define IRTKCU_DEBUG_INTERIM_TIMING(level, section)                          \
     if (debug_time >= level) IRTKCU_INTERIM_TIMING(section)
#else
#  define IRTKCU_DEBUG_INTERIM_TIMING(level, section)   do {} while (false)
#endif

// ---------------------------------------------------------------------------
/// End measurement of execution time of current code block.
///
/// In the following example, the timing measurement are only printed if
/// the global execution time debugging level is greater or equal the
/// specified debugging level. Hence, if debug_time is 1, only the time of
/// the following example section is printed. If debug_time is greater than 1,
/// also the times needed for each execution of inner block are printed together
/// with the summarizing total execution time of the example section. Otherwise,
/// if debug_time is less than 1, no time measurements are printed at all.
/// @code
/// {
///   IRTKCU_START_TIMING();
///   // launch CUDA kernel here
///   IRTKCU_DEBUG_TIMING(1, "kernel name");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at runtime depending on the global
///       variable debug_time.
///
/// @sa IRTKCU_START_TIMING
#ifdef USE_TIMING
#  define IRTKCU_DEBUG_TIMING(level, section)                                  \
     if (debug_time >= level) IRTKCU_END_TIMING(section);                      \
     else do {                                                                 \
       CudaSafeCall( cudaEventDestroy(e_start) );                              \
       CudaSafeCall( cudaEventDestroy(e_stop) );                               \
     } while (false)
#else
#  define IRTKCU_DEBUG_TIMING(level, section)   do {} while (false)
#endif


#endif
