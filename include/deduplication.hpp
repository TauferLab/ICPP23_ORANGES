#ifndef DEDUPLICATION_HPP
#define DEDUPLICATION_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_Bitset.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include "hash_functions.hpp"
#include <fstream>
#include <chrono>
#include <vector>
#include "utils.hpp"

//#define PRINT_CHANGES
#define COALESCE_ON_GPU
#define MEMCPY_GATHER

size_t NUM_CHECKPOINTS=1;

template<typename ExecSpace, typename T>
KOKKOS_INLINE_FUNCTION
bool identical_hashes(const T* a, const T* b, size_t len) {
  if(std::is_same<ExecSpace, Kokkos::DefaultHostExecutionSpace>::value) {
    return memcmp(a, b, len*sizeof(T)) == 0;
  } else {
    for(size_t i=0; i<len; i++) {
      if(a[i] != b[i]) {
        return false;
      }
    }
    return true;
  }
}

// Calculate hashes
template<typename ExecSpace, typename GDVView, typename HashView>
void calculate_hashes(GDVView gdvs, 
                      HashView hashes, 
                      const int chunk_size, 
                      const int num_hashes) {
  auto hash_policy = Kokkos::RangePolicy<ExecSpace>(0, num_hashes);
  Kokkos::parallel_for("Compute hashes", hash_policy, KOKKOS_LAMBDA(const size_t idx) {
    HASH_FUNC hasher;
    int block_size = chunk_size;
    if(chunk_size*(idx+1) > gdvs.span()*sizeof(uint32_t))
      block_size = (gdvs.span()*sizeof(uint32_t))-(idx*chunk_size);
    hasher.hash((uint8_t*)(gdvs.data()) + (idx*chunk_size), block_size, (uint8_t*)(hashes.data())+idx*hasher.digest_size());
  });
}

template<typename ExecSpace, typename GDVView, typename HashView>
void calculate_hashes_omp(GDVView gdvs, 
                      HashView hashes, 
                      const int chunk_size, 
                      const int num_hashes) {
#pragma omp parallel for
  for(size_t idx=0; idx < num_hashes; idx++) {
    HASH_FUNC hasher;
    int block_size = chunk_size;
    if(chunk_size*(idx+1) > gdvs.span()*sizeof(uint32_t))
      block_size = (gdvs.span()*sizeof(uint32_t))-(idx*chunk_size);
    hasher.hash((uint8_t*)(gdvs.data()) + (idx*chunk_size), block_size, (uint8_t*)(hashes.data())+idx*hasher.digest_size());
  }
}

// Find differences by scanning
template<typename ExecSpace, typename GDVView, typename RegionView, typename CounterView>
void find_differences(GDVView previous_gdvs,
                      GDVView current_gdvs, 
                      RegionView changed_regions, 
                      CounterView num_changes, 
                      const int block_len,
                      const int num_blocks) {
  if(std::is_same<ExecSpace,Kokkos::DefaultHostExecutionSpace>::value) {
    auto block_policy = Kokkos::RangePolicy<ExecSpace>(0, num_blocks);
    Kokkos::parallel_for("Scan for changes", block_policy, KOKKOS_LAMBDA(const size_t idx) {
      int block_size = block_len;
      if(block_len*(idx+1) > current_gdvs.span()*sizeof(uint32_t))
        block_size = (current_gdvs.span()*sizeof(uint32_t))-(idx*block_len);
      if(memcmp((uint8_t*)(previous_gdvs.data())+idx*block_len, (uint8_t*)(current_gdvs.data())+idx*block_len, block_size) != 0) {
        size_t n = Kokkos::atomic_fetch_add( &num_changes(0), 1 );
        changed_regions(n) = idx;
      }
    });
  } else {
    auto block_policy = Kokkos::RangePolicy<ExecSpace>(0, num_blocks);
    Kokkos::parallel_for("Scan for changes", block_policy, KOKKOS_LAMBDA(const size_t idx) {
      int block_size = block_len;
      if(block_len*(idx+1) > current_gdvs.span()*sizeof(uint32_t))
        block_size = (current_gdvs.span()*sizeof(uint32_t))-(idx*block_len);
      for(size_t i=0; i<block_size; i++) {
        if(*((uint8_t*)(previous_gdvs.data())+idx*block_len+i) != *((uint8_t*)(previous_gdvs.data())+idx*block_len+i)) {
          size_t n = Kokkos::atomic_fetch_add( &num_changes(0), 1 );
          changed_regions(n) = idx;
          break;
        }
      }
    });
  }
}

// Compare hashes
template<typename ExecSpace, typename HashView, typename RegionView, typename CounterView>
void compare_hashes(HashView prev_hashes, 
                    HashView curr_hashes, 
                    RegionView changed_regions, 
                    CounterView num_changes, 
                    const int hash_len, 
                    const int num_hashes) {
  // Specialize for 64 and 32 bit hash digests
  if(std::is_same<typename HashView::value_type,uint64_t>::value || std::is_same<typename HashView::value_type,uint32_t>::value) {
    auto hash_policy = Kokkos::RangePolicy<ExecSpace>(0, num_hashes);
    Kokkos::parallel_for("Compare hashes", hash_policy, KOKKOS_LAMBDA(const size_t idx) {
      if(prev_hashes(idx) != curr_hashes(idx)) {
        size_t n = Kokkos::atomic_fetch_add( &num_changes(0), 1 );
        changed_regions(n) = idx;
      }
    });
  } else if(std::is_same<typename HashView::value_type, uint8_t>::value) {
    // On host we can use memcmp
    if(std::is_same<ExecSpace,Kokkos::DefaultHostExecutionSpace>::value) {
      auto hash_policy = Kokkos::RangePolicy<ExecSpace>(0, num_hashes);
      Kokkos::parallel_for("Compare hashes", hash_policy, KOKKOS_LAMBDA(const size_t idx) {
        if(memcmp(&prev_hashes(idx*hash_len), &curr_hashes(idx*hash_len), hash_len) != 0) {
          size_t n = Kokkos::atomic_fetch_add( &num_changes(0), 1 );
          changed_regions(n) = idx;
        }
      });
    } else {
    // On device we need to scan byte by byte
      auto hash_policy = Kokkos::RangePolicy<ExecSpace>(0, num_hashes);
      Kokkos::parallel_for("Compare hashes", hash_policy, KOKKOS_LAMBDA(const size_t idx) {
        for(size_t i=0; i<hash_len; i++) {
          if(prev_hashes(idx*hash_len+i) != curr_hashes(idx*hash_len+i)) {
            size_t n = Kokkos::atomic_fetch_add( &num_changes(0), 1 );
            changed_regions(n) = idx;
            break;
          }
        }
      });
    }
  } else {
    printf("ERROR: Hash View type is not uint64_t*, uint32_t*, or uint8_t*\n");
  }
}

template<typename ExecSpace, typename HashView, typename RegionBitset>
void compare_hashes_bitset(HashView prev_hashes, 
                    HashView curr_hashes, 
                    RegionBitset changed_regions, 
                    size_t& num_changes, 
                    const int hash_len, 
                    const int num_hashes) {
//      changed_regions.reset();
//      auto hash_policy = Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic>>(0, num_hashes);
//
//      Kokkos::parallel_for("Find uniques", hash_policy, KOKKOS_LAMBDA(const size_t i) {
//        bool unique = true;
//        size_t dup_idx = 0;
//        for(size_t j=0; j<num_hashes; j++) {
//          if(i != j && identical_hashes<ExecSpace>(&curr_hashes(j*hash_len), &curr_hashes(i*hash_len), hash_len)) {
//            unique = false;
//            dup_idx = j;
//            break;
//          }
//        }
//        if(unique) {
//          changed_regions.set(i);
//        }
//        if(!unique && i<dup_idx) {
//          changed_regions.set(i);
//        }
//      });
//
//      Kokkos::View<size_t*, ExecSpace> unique_chunks("Unique chunks", changed_regions.count());
//      Kokkos::View<size_t[1], ExecSpace> unique_counter("Num unique");
//
//      Kokkos::parallel_for("Get unique indices", hash_policy, KOKKOS_LAMBDA(const size_t i) {
//        if(changed_regions.test(i)) {
//          size_t idx = Kokkos::atomic_fetch_add(&unique_counter(0), 1);
//          unique_chunks(idx) = i;
//        }      
//      });
//
//      Kokkos::RangePolicy<ExecSpace> comp_policy(0, changed_regions.count());
//      Kokkos::parallel_for("Compare hashes", comp_policy, KOKKOS_LAMBDA(const size_t idx) {
//        const size_t i = unique_chunks(idx);
//        for(size_t j=0; j<num_hashes; j++) {
//          if(identical_hashes<ExecSpace>(&prev_hashes(j*hash_len), &curr_hashes(i*hash_len), hash_len)) {
//            changed_regions.reset(i);
//            break;
//          }
//        }
//      });
//      num_changes = changed_regions.count();

//  changed_regions.reset();
std::cout << "Number of hashes: " << num_hashes << std::endl;
Kokkos::View<size_t[1], ExecSpace> comp_counter("Num comparisons");
Kokkos::Experimental::ScatterView<size_t[1]> comp_counter_sv(comp_counter);
Kokkos::View<size_t*> num_searched("num searched", num_hashes);
Kokkos::deep_copy(num_searched, num_hashes);
  using namespace std::chrono;
size_t n_unique = 0;
auto hash_policy = Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic>>(0, num_hashes);
if(changed_regions.count() == 0) {
START_TIMER(find_uniques);
  Kokkos::parallel_for("Find uniques", hash_policy, KOKKOS_LAMBDA(const size_t i) {
auto access = comp_counter_sv.access();
    bool unique = true;
    size_t dup_idx = 0;
    for(size_t j=0; j<num_hashes; j++) {
access(0) += 1;
      if(i != j && identical_hashes<ExecSpace>(&curr_hashes(j*hash_len), &curr_hashes(i*hash_len), hash_len)) {
        unique = false;
        dup_idx = j;
num_searched(i) = j;
        break;
      }
    }
    if(unique) {
      changed_regions.set(i);
    }
    if(!unique && i<dup_idx) {
      changed_regions.set(i);
    }
  });
STOP_TIMER(find_uniques);
std::cout << "Find unique time: " << duration_cast<duration<double>>(find_uniques_end-find_uniques_beg).count() << std::endl;
n_unique = changed_regions.count();
std::cout << "Found " << changed_regions.count() << " unique hashes\n";
//for(int k=0; k<num_searched.size(); k++) {
//  if(changed_regions.test(k)) 
//    std::cout << k << ", " << num_searched(k) << std::endl;
//}
std::cout << "Found " << num_hashes - changed_regions.count() << " duplicate hashes\n";
//for(int k=0; k<num_searched.size(); k++) {
//  if(!changed_regions.test(k)) 
//    std::cout << k << ", " << num_searched(k) << std::endl;
//}
}
//  Kokkos::UnorderedMap<SHA1::digest, size_t, ExecSpace> unique_chunks(num_hashes);
//  Kokkos::parallel_for("Find unique", Kokkos::RangePolicy<ExecSpace>(0, num_hashes), KOKKOS_LAMBDA(const size_t i) {
//    SHA1::digest dig;
//    memcpy(&dig, &curr_hashes(i*hash_len), hash_len);
//    auto result = unique_chunks.insert(dig, i);
//    if(result.success())
//      changed_regions.set(i);
//  });

Kokkos::deep_copy(num_searched, 0);
START_TIMER(comp_hashes);
  Kokkos::parallel_for("Compare hashes", hash_policy, KOKKOS_LAMBDA(const size_t i) {
    if(changed_regions.test(i)) {
auto access = comp_counter_sv.access();
      for(size_t j=0; j<num_hashes; j++) {
access(0) += 1;
        if(identical_hashes<ExecSpace>(&prev_hashes(j*hash_len), &curr_hashes(i*hash_len), hash_len)) {
          changed_regions.reset(i);
num_searched(i) = j;
          break;
        }
      }
    }
  });
STOP_TIMER(comp_hashes);
std::cout << "Compare hashes time: " << duration_cast<duration<double>>(comp_hashes_end-comp_hashes_beg).count() << std::endl;
std::cout << "Found " << n_unique-changed_regions.count() << " old hashes\n";

  num_changes = changed_regions.count();
Kokkos::Experimental::contribute(comp_counter, comp_counter_sv);
printf("Number of comparisons: %u\n", comp_counter(0));
}

// Compare hashes but implemented with OpenMP to prevent blocking the device kernels
template<typename ExecSpace, typename HashView, typename RegionView, typename CounterView>
void compare_hashes_omp(HashView prev_hashes, 
                    HashView curr_hashes, 
                    RegionView changed_regions, 
                    CounterView num_changes, 
                    const int hash_len, 
                    const int num_hashes) {
  if(std::is_same<typename HashView::value_type,uint64_t>::value || std::is_same<typename HashView::value_type,uint32_t>::value) {
#pragma omp parallel for
    for(size_t idx=0; idx<num_hashes; idx++) {
      if(prev_hashes(idx) != curr_hashes(idx)) {
        size_t n = Kokkos::atomic_fetch_add( &num_changes(0), 1 );
        changed_regions(n) = idx;
      }
    }
  } else if(std::is_same<typename HashView::value_type, uint8_t>::value) {
    if(std::is_same<ExecSpace,Kokkos::DefaultHostExecutionSpace>::value) {
#pragma omp parallel for
      for(size_t idx=0; idx<num_hashes; idx++) {
        if(memcmp(&prev_hashes(idx*hash_len), &curr_hashes(idx*hash_len), hash_len) != 0) {
          size_t n = Kokkos::atomic_fetch_add( &num_changes(0), 1 );
          changed_regions(n) = idx;
        }
      }
    } else {
#pragma omp parallel for
      for(size_t idx=0; idx<num_hashes; idx++) {
        for(size_t i=0; i<hash_len; i++) {
          if(prev_hashes(idx*hash_len+i) != curr_hashes(idx*hash_len+i)) {
            size_t n = Kokkos::atomic_fetch_add( &num_changes(0), 1 );
            changed_regions(n) = idx;
            break;
          }
        }
      }
    }
  } else {
    printf("ERROR: Hash View type is not uint64_t*, uint32_t*, or uint8_t*\n");
  }
}

// Gather changes into a contiguous buffer
template <class ExecSpace, typename GDVView, typename RegionView, typename CounterView, typename BufferView>
void gather_changes(GDVView gdvs,
                    RegionView changed_regions,
                    CounterView diff_size,
                    int chunk_size,
                    int num_changes,
                    BufferView diff_buff) {
  auto team_policy = Kokkos::TeamPolicy<ExecSpace>(num_changes, Kokkos::AUTO);
  using team_member_type = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  Kokkos::parallel_for("Gather updates", team_policy, KOKKOS_LAMBDA(team_member_type team_member) {
    size_t idx = team_member.league_rank();
    size_t num_write = chunk_size;
    if(chunk_size*(changed_regions(idx)+1) > gdvs.span()*sizeof(uint32_t))
      num_write = (gdvs.span()*sizeof(uint32_t))-changed_regions(idx)*chunk_size;
    Kokkos::single(Kokkos::PerTeam(team_member), [=] () {
      Kokkos::atomic_add(&diff_size(0), num_write);
    });
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, num_write), [=] (size_t i) {
      *((uint8_t*)(diff_buff.data())+chunk_size*idx+i) = *((uint8_t*)(gdvs.data())+chunk_size*changed_regions(idx)+i);
    });
  });

//  auto diff_policy = Kokkos::RangePolicy<ExecSpace>(0, num_changes);
//  Kokkos::parallel_for("Gather updates", diff_policy, KOKKOS_LAMBDA(const size_t idx) {
//    int num_write = chunk_size;
//    if(chunk_size*(changed_regions(idx)+1) > gdvs.span()*sizeof(uint32_t))
//      num_write = (gdvs.span()*sizeof(uint32_t))-changed_regions(idx)*chunk_size;
//    Kokkos::atomic_add(&diff_size(0), num_write);
//#ifdef MEMCPY_GATHER
//    memcpy((uint8_t*)(diff_buff.data())+chunk_size*idx, (uint8_t*)(gdvs.data())+chunk_size*changed_regions(idx), num_write);
//#else
//    for(size_t i=0; i<num_write; i++) {
//      *((uint8_t*)(diff_buff.data())+chunk_size*idx+i) = *((uint8_t*)(gdvs.data())+chunk_size*changed_regions(idx)+i);
//    }
//#endif
//  });
}

template <class ExecSpace, typename GDVView, typename RegionView, typename CounterView, typename BufferView>
void gather_changes_bitset(GDVView gdvs,
                    RegionView changed_regions,
                    CounterView diff_size,
                    int chunk_size,
                    int num_changes,
                    int num_hashes,
                    BufferView diff_buff) {
  Kokkos::View<size_t[1], ExecSpace> counter("Counter");
  Kokkos::deep_copy(counter, 0);
  Kokkos::View<size_t*, ExecSpace> updated_regions("Changed regions", num_changes);
  Kokkos::parallel_for("Convert bitmap", Kokkos::RangePolicy<ExecSpace>(0, num_hashes), KOKKOS_LAMBDA(const size_t i) {
    if(changed_regions.test(i)) {
      size_t pos = Kokkos::atomic_fetch_add(&counter(0), 1);
      updated_regions(pos) = i;
    }
  });
  Kokkos::fence();
  auto team_policy = Kokkos::TeamPolicy<ExecSpace>(num_changes, Kokkos::AUTO);
  using team_member_type = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  Kokkos::parallel_for("Gather updates", team_policy, KOKKOS_LAMBDA(team_member_type team_member) {
    size_t idx = team_member.league_rank();
    int num_write = chunk_size;
    if(chunk_size*(updated_regions(idx)+1) > gdvs.span()*sizeof(uint32_t))
      num_write = (gdvs.span()*sizeof(uint32_t))-updated_regions(idx)*chunk_size;
    Kokkos::single(Kokkos::PerTeam(team_member), [=] () {
      Kokkos::atomic_add(&diff_size(0), static_cast<size_t>(num_write));
    });
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, num_write), [=] (size_t i) {
      *((uint8_t*)(diff_buff.data())+chunk_size*idx+i) = *((uint8_t*)(gdvs.data())+chunk_size*updated_regions(idx)+i);
    });
  });

//  auto diff_policy = Kokkos::RangePolicy<ExecSpace>(0, num_changes);
//  Kokkos::parallel_for("Gather updates", diff_policy, KOKKOS_LAMBDA(const size_t idx) {
//    int num_write = chunk_size;
//    if(chunk_size*(changed_regions(idx)+1) > gdvs.span()*sizeof(uint32_t))
//      num_write = (gdvs.span()*sizeof(uint32_t))-changed_regions(idx)*chunk_size;
//    Kokkos::atomic_add(&diff_size(0), num_write);
//#ifdef MEMCPY_GATHER
//    memcpy((uint8_t*)(diff_buff.data())+chunk_size*idx, (uint8_t*)(gdvs.data())+chunk_size*changed_regions(idx), num_write);
//#else
//    for(size_t i=0; i<num_write; i++) {
//      *((uint8_t*)(diff_buff.data())+chunk_size*idx+i) = *((uint8_t*)(gdvs.data())+chunk_size*changed_regions(idx)+i);
//    }
//#endif
//  });
}

// Gather changes into a contiguous buffer with OpenMP to avoid blocking device kernels
template <typename ExecSpace, typename GDVView, typename RegionView, typename CounterView, typename BufferView>
void gather_changes_omp(GDVView gdvs,
                    RegionView changed_regions,
                    CounterView diff_size,
                    int chunk_size,
                    int num_changes,
                    BufferView diff_buff) {
#pragma omp parallel for
  for(size_t idx=0; idx<num_changes; idx++) {
    int num_write = chunk_size;
    if(chunk_size*(changed_regions(idx)+1) > gdvs.span()*sizeof(uint32_t))
      num_write = (gdvs.span()*sizeof(uint32_t))-changed_regions(idx)*chunk_size;
    Kokkos::atomic_add(&diff_size(0), num_write);
    memcpy((uint8_t*)(diff_buff.data())+chunk_size*idx, (uint8_t*)(gdvs.data())+chunk_size*changed_regions(idx), num_write);
  }
}

void dedup_case0( GDVs current_gdvs,
                  GDVs::HostMirror previous_gdvs_h,
                  Kokkos::View<size_t*> changed_regions,
                  Kokkos::View<size_t*>::HostMirror changed_regions_h,
                  Kokkos::View<uint8_t*, Kokkos::HostSpace> diff_buff_h,
                  std::string logfile,
                  int chunk_size,
                  int num_blocks,
                  int i
                                                                      ) 
{
  using namespace std::chrono;
  Kokkos::View<size_t[1]> num_changes("Number of changes");
  Kokkos::View<size_t[1]>::HostMirror num_changes_h = Kokkos::create_mirror_view(num_changes);
  Kokkos::deep_copy(num_changes_h, 0);
  Kokkos::deep_copy(changed_regions_h, 0);
  auto current_gdvs_h = Kokkos::create_mirror_view(current_gdvs);
  std::fstream fs, changes;
  fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
#ifdef PRINT_CHANGES
  changes.open(logfile+std::string(".updated_regions"), std::fstream::out | std::fstream::app);
#endif
  Kokkos::fence();
  if(i == 0) {
    START_TIMER(dedup);
    // Copy GDVs to CPU
    START_TIMER(copy);
    Kokkos::deep_copy(previous_gdvs_h, current_gdvs);
    STOP_TIMER(copy);
    STOP_TIMER(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i << ", " 
        << num_blocks << ", "
        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
        << "0" << ", " 
        << "0" << ", " 
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
#ifdef PRINT_CHANGES
    for(size_t i=0; i<num_blocks; i++) {
      changes << i << " ";
    }
    changes << "\n";
#endif
  } else if(i > 0) {
    START_TIMER(dedup);
    // Copy GDVs to CPU
    START_TIMER(copy);
    Kokkos::deep_copy(current_gdvs_h, current_gdvs);
    STOP_TIMER(copy);
    // Find differences
    Kokkos::View<size_t[1]> diff_size("Diff size");
    Kokkos::View<size_t[1]>::HostMirror diff_size_h = Kokkos::create_mirror_view(diff_size);
    START_TIMER(comp);
    find_differences<Kokkos::DefaultHostExecutionSpace>(previous_gdvs_h, current_gdvs_h, changed_regions_h, num_changes_h, chunk_size, num_blocks);
    STOP_TIMER(comp);
    START_TIMER(create_diff);
    gather_changes<Kokkos::DefaultHostExecutionSpace>(current_gdvs_h, changed_regions_h, diff_size_h, chunk_size, num_changes_h(0), diff_buff_h);
    STOP_TIMER(create_diff);
    STOP_TIMER(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i << ", " 
        << num_changes_h(0) << ", "
        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
        << duration_cast<duration<double>>(comp_end-comp_beg).count() << ", " 
        << duration_cast<duration<double>>(create_diff_end-create_diff_beg).count() << ", " 
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
    Kokkos::deep_copy(previous_gdvs_h, current_gdvs_h);
#ifdef PRINT_CHANGES
    Kokkos::sort(changed_regions_h, 0, num_changes_h(0));
    for(size_t i=0; i<num_changes_h(0); i++) {
      changes << changed_regions_h(i) << " ";
    }
    changes << "\n";
#endif
  }
  fs.close();
#ifdef PRINT_CHANGES
  changes.close();
#endif
}

template<typename DeviceHashView, typename HostHashView>
void dedup_case1( 
                  GDVs graph_GDV,
                  GDVs::HostMirror graph_GDV_h,
                  DeviceHashView prev_hashes,
                  HostHashView prev_hashes_h,
                  DeviceHashView curr_hashes,
                  HostHashView curr_hashes_h,
                  Kokkos::View<size_t*> changed_regions,
                  Kokkos::View<size_t*>::HostMirror changed_regions_h,
                  Kokkos::View<uint8_t*, Kokkos::HostSpace> diff_buff_h,
                  std::string logfile,
                  int chunk_size,
                  int num_hashes,
                  int i,
                  Kokkos::View<HASH_FUNC::DIGEST_TYPE**, Kokkos::DefaultHostExecutionSpace> hist_hashes
                                                                      ) 
{
  using namespace std::chrono;
  Kokkos::View<size_t[1]> num_changes("Number of changes");
  Kokkos::View<size_t[1]>::HostMirror num_changes_h = Kokkos::create_mirror_view(num_changes);
  Kokkos::deep_copy(num_changes_h, 0);
  Kokkos::deep_copy(changed_regions_h, 0);
  HASH_FUNC hasher;
  std::fstream fs, changes;
  fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
#ifdef PRINT_CHANGES
  changes.open(logfile+std::string(".updated_regions"), std::fstream::out | std::fstream::app);
#endif
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
    STOP_TIMER(hash);
    Kokkos::deep_copy(prev_hashes_h, curr_hashes_h);
//    auto prev_hashes_subview = Kokkos::subview(hist_hashes, i, Kokkos::ALL);
//    Kokkos::deep_copy(prev_hashes_subview, curr_hashes_h);
    STOP_TIMER(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i << ", " 
        << num_hashes << ", "
        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", "
        << "0" << ", " 
        << "0" << ", " 
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
#ifdef PRINT_CHANGES
    for(size_t i=0; i<num_hashes; i++) {
      changes << i << " ";
    }
    changes << "\n";
#endif
  } else if(i > 0) {
    START_TIMER(dedup);
    // Copy GDVs to CPU
    START_TIMER(copy);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER(copy);
    // Compute hashes on the CPU
    START_TIMER(hash);
    calculate_hashes<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, curr_hashes_h, chunk_size, num_hashes);
    STOP_TIMER(hash);
    // Find differences
    Kokkos::View<size_t[1]> diff_size("Diff size");
    Kokkos::View<size_t[1]>::HostMirror diff_size_h = Kokkos::create_mirror_view(diff_size);
    Kokkos::Bitset<Kokkos::DefaultHostExecutionSpace> updated_regions(num_hashes);
    updated_regions.reset();
    START_TIMER(comp);
    size_t num_updates = 0;
//    uint32_t start = i<NUM_CHECKPOINTS ? 0 : i%NUM_CHECKPOINTS;
//    uint32_t n_iters = i<NUM_CHECKPOINTS ? i : NUM_CHECKPOINTS;
//for(int chkpt=0; chkpt<n_iters; chkpt++) {
//  auto subview = Kokkos::subview(hist_hashes, (start+chkpt)%NUM_CHECKPOINTS, Kokkos::ALL);
//  Kokkos::deep_copy(prev_hashes_h, subview);
    compare_hashes_bitset<Kokkos::DefaultHostExecutionSpace>(prev_hashes_h, curr_hashes_h, updated_regions, num_updates, hasher.digest_size(), num_hashes);
//    compare_hashes_bitset<Kokkos::DefaultHostExecutionSpace>(subview, curr_hashes_h, updated_regions, num_updates, hasher.digest_size(), num_hashes);
//}
    STOP_TIMER(comp);
    START_TIMER(create_diff);
    gather_changes_bitset<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, updated_regions, diff_size_h, chunk_size, num_updates, num_hashes, diff_buff_h);
    STOP_TIMER(create_diff);
    STOP_TIMER(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i << ", " 
        << num_updates << ", "
        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", "
        << duration_cast<duration<double>>(comp_end-comp_beg).count() << ", " 
        << duration_cast<duration<double>>(create_diff_end-create_diff_beg).count() << ", " 
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
//uint32_t replace = i%NUM_CHECKPOINTS;
//auto subview = Kokkos::subview(hist_hashes, replace, Kokkos::ALL);
//Kokkos::deep_copy(subview, curr_hashes_h);
    Kokkos::deep_copy(prev_hashes_h, curr_hashes_h);
#ifdef PRINT_CHANGES
    Kokkos::sort(changed_regions_h, 0, num_updates);
    for(size_t i=0; i<num_updates; i++) {
      changes << changed_regions_h(i) << " ";
    }
    changes << "\n";
#endif
  }
  fs.close();
#ifdef PRINT_CHANGES
  changes.close();
#endif
}

template<typename DeviceHashView, typename HostHashView>
void dedup_case1_async( GDVs graph_GDV,
                  GDVs::HostMirror graph_GDV_h,
                  DeviceHashView prev_hashes,
                  HostHashView prev_hashes_h,
                  DeviceHashView curr_hashes,
                  HostHashView curr_hashes_h,
                  Kokkos::View<size_t*> changed_regions,
                  Kokkos::View<size_t*>::HostMirror changed_regions_h,
                  Kokkos::View<uint8_t*, Kokkos::HostSpace> diff_buff_h,
                  std::string logfile,
                  int chunk_size,
                  int num_hashes,
                  int num_iters,
                  int i
                                                                      ) 
{
  using namespace std::chrono;
  Kokkos::View<size_t[1]> num_changes("Number of changes");
  Kokkos::View<size_t[1]>::HostMirror num_changes_h = Kokkos::create_mirror_view(num_changes);
  Kokkos::View<size_t[1]> diff_size("Diff size");
  Kokkos::View<size_t[1]>::HostMirror diff_size_h = Kokkos::create_mirror_view(diff_size);
  HASH_FUNC hasher;
  std::fstream fs, changes;
  fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
#ifdef PRINT_CHANGES
  changes.open(logfile+std::string(".updated_regions"), std::fstream::out | std::fstream::app);
#endif
  if(i == 0) {
    START_TIMER_ASYNC(dedup);
    Kokkos::fence();
    // Copy GDVs to CPU
    START_TIMER_ASYNC(copy);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER_ASYNC(copy);
    STOP_TIMER_ASYNC(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i << ", " 
        << num_hashes << ", "
        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
        << "0" << ", "
        << "0" << ", " 
        << "0" << ", " 
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
  } else if(i == 1) {
    START_TIMER_ASYNC(dedup);
    // Compute hashes on the CPU
    START_TIMER_ASYNC(hash);
    calculate_hashes_omp<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, curr_hashes_h, chunk_size, num_hashes);
    STOP_TIMER_ASYNC(hash);
    Kokkos::fence();
    // Copy GDVs to CPU
    START_TIMER_ASYNC(copy);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER_ASYNC(copy);
    STOP_TIMER_ASYNC(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i-1 << ", " 
        << num_hashes << ", "
        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", "
        << "0" << ", " 
        << "0" << ", " 
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
//    Kokkos::fence();
  } else if(i == num_iters-1) {
    {
    START_TIMER_ASYNC(dedup);
    // Compute hashes on the CPU
    START_TIMER_ASYNC(hash);
    calculate_hashes_omp<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, curr_hashes_h, chunk_size, num_hashes);
    STOP_TIMER_ASYNC(hash);
    // Find differences
    START_TIMER_ASYNC(comp);
    compare_hashes_omp<Kokkos::DefaultHostExecutionSpace>(prev_hashes_h, curr_hashes_h, changed_regions_h, num_changes_h, hasher.digest_size(), num_hashes);
    STOP_TIMER_ASYNC(comp);
    // Gather changes
    START_TIMER_ASYNC(create_diff);
    gather_changes_omp<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, changed_regions_h, diff_size_h, chunk_size, num_changes_h(0), diff_buff_h);
    STOP_TIMER_ASYNC(create_diff);

    Kokkos::fence();
    // Copy GDVs to CPU
    START_TIMER_ASYNC(copy);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER_ASYNC(copy);
    STOP_TIMER_ASYNC(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i-1 << ", " 
        << num_changes_h(0) << ", "
        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", "
        << duration_cast<duration<double>>(comp_end-comp_beg).count() << ", " 
        << duration_cast<duration<double>>(create_diff_end-create_diff_beg).count() << ", " 
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
//    Kokkos::fence();
#ifdef PRINT_CHANGES
      Kokkos::sort(changed_regions_h, 0, num_changes_h(0));
      for(size_t i=0; i<num_changes_h(0); i++) {
        changes << changed_regions_h(i) << " ";
      }
      changes << "\n";
#endif
    }
    Kokkos::deep_copy(num_changes_h, 0);
    Kokkos::deep_copy(prev_hashes_h, curr_hashes_h);
    {
    START_TIMER_ASYNC(dedup);
    // Compute hashes on the CPU
    START_TIMER_ASYNC(hash);
    calculate_hashes_omp<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, curr_hashes_h, chunk_size, num_hashes);
    STOP_TIMER_ASYNC(hash);
    // Find differences
    START_TIMER_ASYNC(comp);
    compare_hashes_omp<Kokkos::DefaultHostExecutionSpace>(prev_hashes_h, curr_hashes_h, changed_regions_h, num_changes_h, hasher.digest_size(), num_hashes);
    STOP_TIMER_ASYNC(comp);
    // Gather changes
    START_TIMER_ASYNC(create_diff);
    gather_changes_omp<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, changed_regions_h, diff_size_h, chunk_size, num_changes_h(0), diff_buff_h);
    STOP_TIMER_ASYNC(create_diff);

    Kokkos::fence();
    // Copy GDVs to CPU
    START_TIMER_ASYNC(copy);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER_ASYNC(copy);
    STOP_TIMER_ASYNC(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i << ", " 
        << num_changes_h(0) << ", "
        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", "
        << duration_cast<duration<double>>(comp_end-comp_beg).count() << ", " 
        << duration_cast<duration<double>>(create_diff_end-create_diff_beg).count() << ", " 
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
//    Kokkos::fence();
#ifdef PRINT_CHANGES
      Kokkos::sort(changed_regions_h, 0, num_changes_h(0));
      for(size_t i=0; i<num_changes_h(0); i++) {
        changes << changed_regions_h(i) << " ";
      }
      changes << "\n";
#endif
    }
  } else if(i > 0) {
    START_TIMER_ASYNC(dedup);
    // Compute hashes on the CPU
    START_TIMER_ASYNC(hash);
    calculate_hashes_omp<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, curr_hashes_h, chunk_size, num_hashes);
    STOP_TIMER_ASYNC(hash);
    // Find differences
    START_TIMER_ASYNC(comp);
    compare_hashes_omp<Kokkos::DefaultHostExecutionSpace>(prev_hashes_h, curr_hashes_h, changed_regions_h, num_changes_h, hasher.digest_size(), num_hashes);
    STOP_TIMER_ASYNC(comp);
    // Gather changes
    START_TIMER_ASYNC(create_diff);
    gather_changes_omp<Kokkos::DefaultHostExecutionSpace>(graph_GDV_h, changed_regions_h, diff_size_h, chunk_size, num_changes_h(0), diff_buff_h);
    STOP_TIMER_ASYNC(create_diff);

    Kokkos::fence();
    // Copy GDVs to CPU
    START_TIMER_ASYNC(copy);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER_ASYNC(copy);
    STOP_TIMER_ASYNC(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i << ", " 
        << num_changes_h(0) << ", "
        << duration_cast<duration<double>>(copy_end-copy_beg).count() << ", " 
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", "
        << duration_cast<duration<double>>(comp_end-comp_beg).count() << ", " 
        << duration_cast<duration<double>>(create_diff_end-create_diff_beg).count() << ", " 
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
//    Kokkos::fence();
#endif
#ifdef PRINT_CHANGES
    Kokkos::sort(changed_regions_h, 0, num_changes_h(0));
    for(size_t i=0; i<num_changes_h(0); i++) {
      changes << changed_regions_h(i) << " ";
    }
    changes << "\n";
#endif
  }
  fs.close();
#ifdef PRINT_CHANGES
  changes.close();
#endif
  Kokkos::deep_copy(prev_hashes_h, curr_hashes_h);
}

template<typename DeviceHashView, typename HostHashView>
void dedup_case2( 
                  GDVs graph_GDV,
                  GDVs::HostMirror graph_GDV_h,
                  DeviceHashView prev_hashes,
                  HostHashView prev_hashes_h,
                  DeviceHashView curr_hashes,
                  HostHashView curr_hashes_h,
                  Kokkos::View<size_t*> changed_regions,
                  Kokkos::View<size_t*>::HostMirror changed_regions_h,
                  Kokkos::View<uint8_t*, Kokkos::HostSpace> diff_buff_h,
                  std::string logfile,
                  int chunk_size,
                  int num_hashes,
                  int i,
                  Kokkos::View<HASH_FUNC::DIGEST_TYPE**, Kokkos::DefaultHostExecutionSpace> hist_hashes
                                                                      ) 
{
  using namespace std::chrono;
  Kokkos::View<size_t[1]> num_changes("Number of changes");
  Kokkos::View<size_t[1]>::HostMirror num_changes_h = Kokkos::create_mirror_view(num_changes);
  Kokkos::deep_copy(num_changes_h, 0);
  Kokkos::fence();
  HASH_FUNC hasher;
  Kokkos::deep_copy(changed_regions_h, 0);
  std::fstream fs, changes;
  fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
#ifdef PRINT_CHANGES
  changes.open(logfile+std::string(".updated_regions"), std::fstream::out | std::fstream::app);
#endif
  Kokkos::fence();
  if(i == 0) {
    START_TIMER(dedup);
    // Compute hashes on the GPU
    START_TIMER(hash);
    calculate_hashes<Kokkos::DefaultExecutionSpace>(graph_GDV, curr_hashes, chunk_size, num_hashes);
    STOP_TIMER(hash);
    START_TIMER(copy_diff);
    Kokkos::deep_copy(prev_hashes_h, curr_hashes);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER(copy_diff);
    STOP_TIMER(dedup);
#ifdef PRINT_TIMERS
    fs  << i << ", " 
        << num_hashes << ", "
        << "0" << ", " 
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", " 
        << "0" << ", " 
        << "0" << ", "
        << "0" << ", " 
        << duration_cast<duration<double>>(copy_diff_end-copy_diff_beg).count() << ", "
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
//    Kokkos::deep_copy(prev_hashes_h, curr_hashes);
    auto prev_hashes_subview = Kokkos::subview(hist_hashes, i, Kokkos::ALL);
    Kokkos::deep_copy(prev_hashes_subview, curr_hashes);
#ifdef PRINT_CHANGES
    for(size_t i=0; i<num_hashes; i++) {
      changes << i << " ";
    }
    changes << "\n";
#endif
  } else if(i > 0) {
    START_TIMER(dedup);
    // Copy previous hashes to GPU
    START_TIMER(copy_hash);
    Kokkos::deep_copy(prev_hashes, prev_hashes_h);
    STOP_TIMER(copy_hash);
    // Compute hashes on the GPU
    START_TIMER(hash);
    calculate_hashes<Kokkos::DefaultExecutionSpace>(graph_GDV, curr_hashes, chunk_size, num_hashes);
    STOP_TIMER(hash);
    // Find differences and make diff
    Kokkos::Bitset<Kokkos::DefaultExecutionSpace> updated_regions(num_hashes);
    size_t num_updates = 0;
    START_TIMER(comp);
    uint32_t start = i<NUM_CHECKPOINTS ? 0 : i%NUM_CHECKPOINTS;
    uint32_t n_iters = i<NUM_CHECKPOINTS ? i : NUM_CHECKPOINTS;
for(int chkpt=0; chkpt<n_iters; chkpt++) {
  auto subview = Kokkos::subview(hist_hashes, (start+chkpt)%NUM_CHECKPOINTS, Kokkos::ALL);
  Kokkos::deep_copy(prev_hashes, subview);
    compare_hashes_bitset<Kokkos::DefaultExecutionSpace>(prev_hashes, curr_hashes, updated_regions, num_updates, hasher.digest_size(), num_hashes);
}
    num_changes_h(0) = num_updates;
    STOP_TIMER(comp);
//    START_TIMER(comp);
//    compare_hashes<Kokkos::DefaultExecutionSpace>(prev_hashes, curr_hashes, changed_regions, num_changes, hasher.digest_size(), num_hashes);
//    STOP_TIMER(comp);
    START_TIMER(copy_changes);
    Kokkos::deep_copy(changed_regions_h, changed_regions);
    STOP_TIMER(copy_changes);
    // Create diff
//    START_TIMER(create_diff);
    Kokkos::View<uint8_t*> buffer("Buffer", num_changes_h(0)*chunk_size);
    Kokkos::View<size_t[1]> diff_size("Diff size");
    Kokkos::View<size_t[1]>::HostMirror diff_size_h = Kokkos::create_mirror_view(diff_size);
    START_TIMER(create_diff);
    gather_changes_bitset<Kokkos::DefaultExecutionSpace>(graph_GDV, updated_regions, diff_size, chunk_size, num_changes_h(0), num_hashes, buffer);
    Kokkos::deep_copy(diff_size_h, diff_size);
    STOP_TIMER(create_diff);
    // Copy diff
    START_TIMER(copy_diff);
    Kokkos::deep_copy(prev_hashes_h, curr_hashes);
    auto host = Kokkos::subview(diff_buff_h, Kokkos::make_pair((size_t)(0), (diff_size_h(0))));
    auto device = Kokkos::subview(buffer, Kokkos::make_pair((size_t)(0), (diff_size_h(0))));
    Kokkos::deep_copy(host, device);
    STOP_TIMER(copy_diff);
    STOP_TIMER(dedup);
#ifdef PRINT_TIMERS
    // Print stats
    fs  << i << ", " 
        << num_changes_h(0) << ", "
        << duration_cast<duration<double>>(copy_hash_end-copy_hash_beg).count() << ", " 
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", " 
        << duration_cast<duration<double>>(comp_end-comp_beg).count() << ", " 
        << duration_cast<duration<double>>(copy_changes_end-copy_changes_beg).count() << ", "
        << duration_cast<duration<double>>(create_diff_end-create_diff_beg).count() << ", " 
        << duration_cast<duration<double>>(copy_diff_end-copy_diff_beg).count() << ", "
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
//    Kokkos::deep_copy(prev_hashes_h, curr_hashes);
uint32_t replace = i%NUM_CHECKPOINTS;
auto subview = Kokkos::subview(hist_hashes, replace, Kokkos::ALL);
Kokkos::deep_copy(subview, curr_hashes);
#ifdef PRINT_CHANGES
    Kokkos::sort(changed_regions_h, 0, num_changes_h(0));
    for(size_t i=0; i<num_changes_h(0); i++) {
      changes << changed_regions_h(i) << " ";
    }
    changes << "\n";
#endif
  }
  fs.close();
#ifdef PRINT_CHANGES
  changes.close();
#endif
}

template<typename DeviceHashView, typename HostHashView>
void dedup_case3( GDVs graph_GDV,
                  GDVs::HostMirror graph_GDV_h,
                  DeviceHashView prev_hashes,
                  HostHashView prev_hashes_h,
                  DeviceHashView curr_hashes,
                  HostHashView curr_hashes_h,
                  Kokkos::View<size_t*> changed_regions,
                  Kokkos::View<size_t*>::HostMirror changed_regions_h,
                  Kokkos::View<uint8_t*, Kokkos::HostSpace> diff_buff_h,
                  std::string logfile,
                  int chunk_size,
                  int num_hashes,
                  int i
                                                                      ) 
{
  using namespace std::chrono;
  Kokkos::fence();
  Kokkos::View<size_t[1]> num_changes("Number of changes");
  Kokkos::View<size_t[1]>::HostMirror num_changes_h = Kokkos::create_mirror_view(num_changes);
  HASH_FUNC hasher;
  Kokkos::deep_copy(num_changes_h, 0);
  Kokkos::deep_copy(changed_regions_h, 0);
  std::fstream fs, changes;
  fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
#ifdef PRINT_CHANGES
  changes.open(logfile+std::string(".updated_regions"), std::fstream::out | std::fstream::app);
#endif
  Kokkos::fence();
  if(i == 0) {
    START_TIMER(dedup);
    // Compute hashes on the GPU
    START_TIMER(hash);
    calculate_hashes<Kokkos::DefaultExecutionSpace>(graph_GDV, curr_hashes, chunk_size, num_hashes);
    STOP_TIMER(hash);
    // Transfer hashes to CPU
    START_TIMER(copy_hash);
    Kokkos::deep_copy(prev_hashes_h, curr_hashes);
    STOP_TIMER(copy_hash);
    // Copy GDVs to CPU
    START_TIMER(copy_diff);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    STOP_TIMER(copy_diff);
    STOP_TIMER(dedup);
#ifdef PRINT_TIMERS
    fs  << i << ", "  // Chunk
        << num_hashes << ", "
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", " // Compute hashes
        << duration_cast<duration<double>>(copy_hash_end-copy_hash_beg).count() << ", " // Copy hashes
        << "0, " // Compare hashes
        << "0, " // Create diff
        << duration_cast<duration<double>>(copy_diff_end-copy_diff_beg).count() << ", " // Copy diff
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
    Kokkos::deep_copy(prev_hashes_h, curr_hashes);
#ifdef PRINT_CHANGES
    for(size_t i=0; i<num_hashes; i++) {
      changes << i << " ";
    }
    changes << "\n";
#endif
  } else if(i > 0) {
    START_TIMER(dedup);
    // Compute hashes on the GPU
    START_TIMER(hash);
    calculate_hashes<Kokkos::DefaultExecutionSpace>(graph_GDV, curr_hashes, chunk_size, num_hashes);
    STOP_TIMER(hash);
    // Transfer hashes to CPU
    START_TIMER(copy_hash);
    Kokkos::deep_copy(curr_hashes_h, curr_hashes);
    STOP_TIMER(copy_hash);
    // Compare hashes on CPU
    START_TIMER(comp);
    compare_hashes<Kokkos::DefaultHostExecutionSpace>(prev_hashes_h, curr_hashes_h, changed_regions_h, num_changes_h, hasher.digest_size(), num_hashes);
    STOP_TIMER(comp);
    Kokkos::View<uint8_t*> buffer("Buffer", num_changes_h(0)*chunk_size);
    Kokkos::View<size_t[1]> diff_size("Diff size");
    Kokkos::View<size_t[1]>::HostMirror diff_size_h = Kokkos::create_mirror_view(diff_size);
    // Create diff on GPU
    START_TIMER(create_diff);
//    Kokkos::deep_copy(changed_regions, changed_regions_h);
//    gather_changes<Kokkos::DefaultExecutionSpace>(graph_GDV, changed_regions, diff_size, chunk_size, num_changes_h(0), buffer);
//  //  for(int idx=0; idx<num_changes_h(0); idx++) {
//    Kokkos::parallel_for("Create diff", Kokkos::RangePolicy<>(0,num_changes_h(0)), KOKKOS_LAMBDA(const size_t idx) {
//      int num_write = chunk_size;
//  //    if(chunk_size*(changed_regions_h(idx)+1) > graph_GDV.span()*sizeof(uint32_t))
//  //      num_write = graph_GDV.span()*sizeof(uint32_t) - changed_regions_h(idx)*chunk_size;
//  //    diff_size_h(0) += num_write;
//      if(chunk_size*(changed_regions(idx)+1) > graph_GDV.span()*sizeof(uint32_t))
//        num_write = graph_GDV.span()*sizeof(uint32_t) - changed_regions(idx)*chunk_size;
//  //    diff_size(0) += num_write;
//      Kokkos::atomic_add(&diff_size(0), num_write);
//      for(size_t i=0; i<num_write; i++) {
//        buffer(chunk_size*idx+i) = *((uint8_t*)(graph_GDV.data())+chunk_size*changed_regions(idx)+i);
//      }
//  //    cudaMemcpy( (uint8_t*)(buffer.data())+chunk_size*idx, 
//  //                (uint8_t*)(graph_GDV.data())+chunk_size*changed_regions_h(idx), 
//  //                num_write, cudaMemcpyDeviceToDevice);
//    });
//    Kokkos::deep_copy(diff_size_h, diff_size);
  //  }
#ifdef COALESCE_ON_GPU
    gather_changes<Kokkos::DefaultExecutionSpace>(graph_GDV, changed_regions, diff_size, chunk_size, num_changes_h(0), buffer);
#else
    pull_scattered_changes<Kokkos::DefaultHostExecutionSpace>(graph_GDV, changed_regions_h, diff_size_h, chunk_size, num_changes_h(0), diff_buff_h);
#endif
    STOP_TIMER(create_diff);
    // Copy diff to CPU
    START_TIMER(copy_diff);
    auto host = Kokkos::subview(diff_buff_h, Kokkos::make_pair((size_t)(0), (diff_size_h(0))));
    auto device = Kokkos::subview(buffer, Kokkos::make_pair((size_t)(0), (diff_size_h(0))));
    Kokkos::deep_copy(host, device);
    STOP_TIMER(copy_diff);
    STOP_TIMER(dedup);
#ifdef PRINT_TIMERS
    fs  << i << ", " 
        << num_changes_h(0) << ", "
        << duration_cast<duration<double>>(hash_end-hash_beg).count() << ", " 
        << duration_cast<duration<double>>(copy_hash_end-copy_hash_beg).count() << ", " 
        << duration_cast<duration<double>>(comp_end-comp_beg).count() << ", " 
        << duration_cast<duration<double>>(create_diff_end-create_diff_beg).count() << ", " 
        << duration_cast<duration<double>>(copy_diff_end-copy_diff_beg).count() << ", "
        << duration_cast<duration<double>>(dedup_end-dedup_beg).count() << std::endl; // Total time
#endif
    Kokkos::deep_copy(prev_hashes_h, curr_hashes);
#ifdef PRINT_CHANGES
    Kokkos::sort(changed_regions_h, 0, num_changes_h(0));
    for(size_t i=0; i<num_changes_h(0); i++) {
      changes << changed_regions_h(i) << " ";
    }
    changes << "\n";
#endif
  }
  fs.close();
#ifdef PRINT_CHANGES
  changes.close();
#endif
}
#endif
