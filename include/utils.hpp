#ifndef UTILS_HPP
#define UTILS_HPP

#include <Kokkos_Core.hpp>

//#define xstr(s) str(s)
//#define str(s) #s

//#define START_TIMER(timer) Kokkos::Profiling::pushRegion( #timer );  \
//                          Kokkos::fence()
//#define STOP_TIMER(timer) Kokkos::fence(); \
//                          Kokkos::Profiling::popRegion()

#define START_TIMER(timer) auto timer##_beg = high_resolution_clock::now(); \
                          Kokkos::fence();
#define STOP_TIMER(timer) Kokkos::fence();                                  \
                          auto timer##_end = high_resolution_clock::now();

#endif
