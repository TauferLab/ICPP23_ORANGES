#ifndef DEDUPLICATION_HPP
#define DEDUPLICATION_HPP

#include <Kokkos_Core.hpp>
#include "hash_functions.hpp"
#include <fstream>
#include <chrono>
#include "utils.hpp"

// Calculate hashes
template<typename ExecSpace, typename GDVView, typename HashView>
void calculate_hashes(GDVView gdvs, 
                      HashView hashes, 
                      const int chunk_size, 
                      const int num_hashes) {
  auto hash_policy = Kokkos::RangePolicy<ExecSpace>(0, num_hashes);
  Kokkos::parallel_for("Hashes on CPU", hash_policy, KOKKOS_LAMBDA(const size_t idx) {
    int hash_size = chunk_size;
    if(chunk_size*(idx+1) > gdvs.span()*sizeof(uint32_t))
      hash_size = (gdvs.span()*sizeof(uint32_t))-(idx*chunk_size);
    uint32_t hash = hash32((uint8_t*)(gdvs.data())+(idx*chunk_size), hash_size);
    hashes(idx) = hash;
  });
}

// Compare hashes
template<typename ExecSpace, typename HashView, typename RegionView, typename CounterView>
void compare_hashes(HashView prev_hashes, 
                    HashView curr_hashes, 
                    RegionView changed_regions, 
                    CounterView num_changes, 
                    const int chunk_size, 
                    const int num_hashes) {
  auto hash_policy = Kokkos::RangePolicy<ExecSpace>(0, num_hashes);
  Kokkos::parallel_for("Compare hashes on CPU", hash_policy, KOKKOS_LAMBDA(const size_t idx) {
    if(prev_hashes(idx) != curr_hashes(idx)) {
      size_t n = Kokkos::atomic_fetch_add( &num_changes(0), 1 );
      changed_regions(n) = idx;
    }
  });
}

// Gather changes
template <typename ExecSpace, typename GDVView, typename RegionView, typename CounterView, typename BufferView>
void gather_changes(GDVView gdvs,
                    RegionView changed_regions,
                    CounterView diff_size,
                    int chunk_size,
                    int num_changes,
                    BufferView diff_buff) {
  auto diff_policy = Kokkos::RangePolicy<ExecSpace>(0, num_changes);
  Kokkos::parallel_for("Create diff on CPU", diff_policy, KOKKOS_LAMBDA(const size_t idx) {
    int num_write = chunk_size;
    if(chunk_size*(changed_regions(idx)+1) > gdvs.span()*sizeof(uint32_t))
      num_write = (gdvs.span()*sizeof(uint32_t))-changed_regions(idx)*chunk_size;
    Kokkos::atomic_add(&diff_size(0), num_write);
    std::memcpy((uint8_t*)(diff_buff.data())+chunk_size*idx, (uint8_t*)(gdvs.data())+chunk_size*changed_regions(idx), num_write);
  });
}

//template<typename DeviceUint32View, typename HostUint32View, typename DeviceSizetView, typename HostSizetView>
void dedup_case1( GDVs graph_GDV,
                  GDVs::HostMirror graph_GDV_h,
                  Kokkos::View<uint32_t*> prev_hashes,
                  Kokkos::View<uint32_t*>::HostMirror prev_hashes_h,
                  Kokkos::View<uint32_t*> curr_hashes,
                  Kokkos::View<uint32_t*>::HostMirror curr_hashes_h,
                  Kokkos::View<size_t*> changed_regions,
                  Kokkos::View<size_t*>::HostMirror changed_regions_h,
                  Kokkos::View<uint8_t*, Kokkos::HostSpace> diff_buff_h,
                  std::fstream& fs,
                  int chunk_size,
                  int num_hashes,
                  int i
                                                                      ) 
{
  using namespace std::chrono;
  Kokkos::View<size_t[1]> num_changes("Number of changes");
  Kokkos::View<size_t[1]>::HostMirror num_changes_h = Kokkos::create_mirror_view(num_changes);
  Kokkos::deep_copy(num_changes_h, 0);
  Kokkos::deep_copy(changed_regions_h, 0);
  Kokkos::fence();
  if(i == 0) {
    START_TIMER(dedup);
    // Copy GDVs to CPU
    START_TIMER(copy);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER(copy);
    // Compute hashes on the CPU
    START_TIMER(hash);
    calculate_hashes<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, curr_hashes_h, chunk_size, num_hashes);
//    Kokkos::parallel_for("Hashes on CPU", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0,num_hashes), KOKKOS_LAMBDA(const size_t idx) {
//      int hash_size = chunk_size;
//      if(chunk_size*(idx+1) > graph_GDV_h.span()*sizeof(uint32_t))
//        hash_size = (graph_GDV_h.span()*sizeof(uint32_t))-(idx*chunk_size);
//      uint32_t hash = hash32((uint8_t*)(graph_GDV_h.data())+(idx*chunk_size), hash_size);
//      curr_hashes_h(idx) = hash;
//    });
    STOP_TIMER(hash);
    STOP_TIMER(dedup);
    // Print stats
//    fs  << i << ", " 
//        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
//        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", "
//        << "0" << ", " 
//        << "0" << ", " 
//        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
    Kokkos::deep_copy(prev_hashes_h, curr_hashes_h);
  } else if(i > 0) {
    START_TIMER(dedup);
    // Copy GDVs to CPU
    START_TIMER(copy);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER(copy);
    // Compute hashes on the CPU
    START_TIMER(hash);
    calculate_hashes<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, curr_hashes_h, chunk_size, num_hashes);
//    Kokkos::parallel_for("Hashes on CPU", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0,num_hashes), KOKKOS_LAMBDA(const size_t idx) {
//      int hash_size = chunk_size;
//      if(chunk_size*(idx+1) > graph_GDV_h.span()*sizeof(uint32_t))
//        hash_size = (graph_GDV_h.span()*sizeof(uint32_t))-(idx*chunk_size);
//      uint32_t hash = hash32((uint8_t*)(graph_GDV_h.data())+(idx*chunk_size), hash_size);
//      curr_hashes_h(idx) = hash;
//    });
    STOP_TIMER(hash);
    // Find differences
    Kokkos::View<size_t[1]> diff_size("Diff size");
    Kokkos::View<size_t[1]>::HostMirror diff_size_h = Kokkos::create_mirror_view(diff_size);
    START_TIMER(comp);
    compare_hashes<Kokkos::DefaultHostExecutionSpace>(prev_hashes_h, curr_hashes_h, changed_regions_h, num_changes_h, chunk_size, num_hashes);
//    Kokkos::parallel_for("Compare hashes on CPU", Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, num_hashes), KOKKOS_LAMBDA(const size_t idx) {
//      if(prev_hashes_h(idx) != curr_hashes_h(idx)) {
//        size_t n = Kokkos::atomic_fetch_add( &num_changes_h(0), 1 );
//        changed_regions_h(n) = idx;
//      }
//    });
    STOP_TIMER(comp);
    START_TIMER(create_diff);
    gather_changes<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, changed_regions_h, diff_size_h, chunk_size, num_changes_h(0), diff_buff_h);
//    auto diff_policy = Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, num_changes_h(0));
//    Kokkos::parallel_for("Create diff on CPU", diff_policy, KOKKOS_LAMBDA(const size_t idx) {
//      int num_write = chunk_size;
//      if(chunk_size*(changed_regions_h(idx)+1) > graph_GDV_h.span()*sizeof(uint32_t))
//        num_write = (graph_GDV_h.span()*sizeof(uint32_t))-changed_regions_h(idx)*chunk_size;
//      Kokkos::atomic_add(&diff_size_h(0), num_write);
//      std::memcpy((uint8_t*)(diff_buff_h.data())+chunk_size*idx, (uint8_t*)(graph_GDV_h.data())+chunk_size*changed_regions_h(idx), num_write);
//    });
    STOP_TIMER(create_diff);
    STOP_TIMER(dedup);
    // Print stats
//    fs  << i << ", " 
//        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
//        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", "
//        << duration_cast<duration<double>>(comp_end-comp_beg).count() << ", " 
//        << duration_cast<duration<double>>(create_diff_end-create_diff_beg).count() << ", " 
//        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
    Kokkos::deep_copy(prev_hashes_h, curr_hashes_h);
  }
}

#endif
