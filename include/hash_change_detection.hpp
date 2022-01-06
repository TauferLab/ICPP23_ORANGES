#include <chrono>
//#include "etl/crc32.h"
//#include "etl/crc64_ecma.h"
//#include "etl/fnv_1.h"
//#include "etl/jenkins.h"
//#include "etl/pearson.h"
#include "etl/murmur3.h"

template<class GraphletDegreeVector>
KOKKOS_INLINE_FUNCTION
void generate_hashes(GraphletDegreeVector gdv, Kokkos::View<uint64_t*>& hashes, size_t blocksize) {
  size_t blocksize_elements = blocksize/4;
  if(hashes.size() != gdv.span()/(blocksize_elements)) {
    printf("Hashes destination does not match number of blocks\n");
  }

  etl::murmur3<uint64_t> murmur3_64_gen;

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//  for(int i=0; i<gdv.span()/(blocksize_elements); i++) {
////    murmur3_64_gen.add((uint8_t*)(gdv.data()+(i*blocksize_elements)), (uint8_t*)(gdv.data()+((i+1)*blocksize_elements)));
//    murmur3_64_gen.add((uint8_t*)(&(gdv(i*blocksize_elements))), (uint8_t*)(&(gdv((i+1)*blocksize_elements))));
//    hashes(i) = murmur3_64_gen.value();
//  }
  Kokkos::parallel_for("Generate hashes", Kokkos::RangePolicy<>(0,gdv.span()/blocksize_elements), KOKKOS_LAMBDA(const int i) {
    etl::murmur3<uint64_t> murmur3_64_gen;
    murmur3_64_gen.add((uint8_t*)(gdv.data()+(i*blocksize_elements)), (uint8_t*)(gdv.data()+((i+1)*blocksize_elements)));
    hashes(i) = murmur3_64_gen.value();
  });
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
//  printf("murmur3_64 hash (blocksize=%d) took \t\t%lf seconds\n", blocksize, time_span.count());
}
