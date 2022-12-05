#ifndef UTILS_HPP
#define UTILS_HPP

#include <Kokkos_Core.hpp>

#define PRINT_TIMERS

#ifdef PRINT_TIMERS
  #define START_TIMER(timer) auto timer##_beg = high_resolution_clock::now(); \
                            Kokkos::fence();
  #define STOP_TIMER(timer) Kokkos::fence();                                  \
                            auto timer##_end = high_resolution_clock::now();
  
  #define START_TIMER_ASYNC(timer) auto timer##_beg = high_resolution_clock::now(); 
  #define STOP_TIMER_ASYNC(timer)  auto timer##_end = high_resolution_clock::now();
#else
  #define START_TIMER(timer) Kokkos::Profiling::pushRegion( #timer );  \
                            Kokkos::fence()
  #define STOP_TIMER(timer) Kokkos::fence(); \
                            Kokkos::Profiling::popRegion()
#endif

#endif
