#pragma once

#include <chrono>
#include <iostream>

// [[nodiscard]] attributes on STL functions
#ifndef _HAS_NODISCARD
#ifndef __has_cpp_attribute
#define _HAS_NODISCARD 0
#elif __has_cpp_attribute(nodiscard) >= \
    201603L  // TRANSITION, VSO#939899 (need toolset update)
#define _HAS_NODISCARD 1
#else
#define _HAS_NODISCARD 0
#endif
#endif  // _HAS_NODISCARD

#if _HAS_NODISCARD
#define _NODISCARD [[nodiscard]]
#else  // ^^^ CAN HAZ [[nodiscard]] / NO CAN HAZ [[nodiscard]] vvv
#define _NODISCARD
#endif  // _HAS_NODISCARD

template <typename Func>
void TimeTask(const std::string& task_name, Func&& f) {
  const auto t0 = std::chrono::high_resolution_clock::now();

  std::forward<Func>(f)();

  const auto t1 = std::chrono::high_resolution_clock::now();
  const auto time_span =
      std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0);

  std::cout << "Finished " << task_name << "! Time took: " << time_span.count()
            << "s. " << std::endl;
}

#define BASE_BITS 8
#define BASE (1 << BASE_BITS)
#define MASK (BASE-1)
#define DIGITS(v, shift) (((v) >> shift) & MASK)
 
void omp_lsd_radix_sort(size_t n, std::vector<uint64_t>& data) {
    std::vector<uint64_t> buffer(n);
    int total_digits = sizeof(uint64_t)*8;
 
    //Each thread use local_bucket to move data
    size_t i;
    for(int shift = 0; shift < total_digits; shift+=BASE_BITS) {
        size_t bucket[BASE] = {0};
 
        size_t local_bucket[BASE] = {0}; // size needed in each bucket/thread
        //1st pass, scan whole and check the count
        #pragma omp parallel firstprivate(local_bucket)
        {
            #pragma omp for schedule(static) nowait
            for(i = 0; i < n; i++){
                local_bucket[DIGITS(data[i], shift)]++;
            }
            #pragma omp critical
            for(i = 0; i < BASE; i++) {
                bucket[i] += local_bucket[i];
            }
            #pragma omp barrier
            #pragma omp single
            for (i = 1; i < BASE; i++) {
                bucket[i] += bucket[i - 1];
            }
            int nthreads = omp_get_num_threads();
            int tid = omp_get_thread_num();
            for(int cur_t = nthreads - 1; cur_t >= 0; cur_t--) {
                if(cur_t == tid) {
                    for(i = 0; i < BASE; i++) {
                        bucket[i] -= local_bucket[i];
                        local_bucket[i] = bucket[i];
                    }
                } else { //just do barrier
                    #pragma omp barrier
                }
 
            }
            #pragma omp for schedule(static)
            for(i = 0; i < n; i++) { //note here the end condition
                buffer[local_bucket[DIGITS(data[i], shift)]++] = data[i];
            }
        }
        //now move data
        std::vector<uint64_t> tmp = data;
        data = buffer;
        buffer = tmp;
    }
}