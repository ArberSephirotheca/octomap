#pragma once

#include <chrono>
#include <iostream>
//#include <librealsense2/rs.hpp>
#include <pcl/point_types.h>
#include <pcl/filters/passthrough.h>

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
float TimeTask(const std::string& task_name, Func&& f) {
  const auto t0 = std::chrono::high_resolution_clock::now();

  std::forward<Func>(f)();

  const auto t1 = std::chrono::high_resolution_clock::now();
  const auto time_span =
      std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0);

  std::cout << "Finished " << task_name << "! Time took: " << time_span.count()
            << "s. " << std::endl;
  return time_span.count();
}

#define BASE_BITS 8
#define BASE (1 << BASE_BITS)
#define MASK (BASE-1)
#define DIGITS(v, shift) (((v) >> shift) & MASK)
 
void omp_lsd_radix_sort(size_t n, std::vector<uint64_t>& data, const int num_threads) {
    std::vector<uint64_t> buffer(n);
    int total_digits = sizeof(uint64_t)*8;
 
    //Each thread use local_bucket to move data
    size_t i;
    for(int shift = 0; shift < total_digits; shift+=BASE_BITS) {
        size_t bucket[BASE] = {0};
 
        size_t local_bucket[BASE] = {0}; // size needed in each bucket/thread
         omp_set_num_threads(num_threads);
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


template <uint8_t Axis>
bool CompareAxis(const Eigen::Vector3f& a, const Eigen::Vector3f& b) {
  if constexpr (Axis == 0) {
    return a.x() < b.x();
  } else if constexpr (Axis == 1) {
    return a.y() < b.y();
  } else {
    return a.z() < b.z();
  }
}

template <uint8_t Axis>
bool CompareAxis_pcl(const pcl::PointXYZ& a, const pcl::PointXYZ& b) {
  if constexpr (Axis == 0) {
    return a.x < b.x;
  } else if constexpr (Axis == 1) {
    return a.y < b.y;
  } else {
    return a.z < b.z;
  }
}
// Function to perform counting sort on a specific digit's place (0 or 1)
void countingSort(std::vector<Code_t>& arr, int exp, int threadId, int numThreads) {
    const size_t n = arr.size();
    std::vector<Code_t> output(n);
    std::vector<Code_t> count(2, 0);

    for (size_t i = threadId; i < n; i += numThreads) {
        count[(arr[i] >> exp) & 1]++;
    }

    // Synchronize threads before updating the count array
    std::mutex mtx;
    std::unique_lock<std::mutex> lock(mtx);
    for (int i = 0; i < numThreads; i++) {
        count[0] += count[2 * i];
        count[1] += count[2 * i + 1];
    }
    lock.unlock();

    for (size_t i = threadId; i < n; i += numThreads) {
        output[count[(arr[i] >> exp) & 1]++] = arr[i];
    }

    for (size_t i = threadId; i < n; i += numThreads) {
        arr[i] = output[i];
    }
}

// Function to perform parallel radix sort using multiple threads
void parallelRadixSort(std::vector<Code_t>& arr, int numThreads) {
    //const size_t n = arr.size();
    const int numBits = 64; // Assuming 64-bit Morton codes

    for (int exp = 0; exp < numBits; exp++) {
        std::vector<std::thread> threads(numThreads);

        for (int i = 0; i < numThreads; i++) {
            threads[i] = std::thread(countingSort, std::ref(arr), exp, i, numThreads);
        }

        for (int i = 0; i < numThreads; i++) {
            threads[i].join();
        }
    }
}

void compute_morton_code_threaded(const int input_size, std::vector<Eigen::Vector3f>& inputs,std::vector<Code_t>& morton_keys, float min_coord, const float range, int num_threads){

	const auto elements_per_thread = math::divide_ceil<int>(input_size, num_threads);
  	const auto worker_fn = [&inputs, &morton_keys,input_size, min_coord, range, elements_per_thread](int i) {
      for(int t = i * elements_per_thread; t < math::min(input_size, (i + 1)*elements_per_thread); ++t){
        Eigen::Vector3f input = inputs[t];
        morton_keys[t] = PointToCode(input.x(), input.y(), input.z(), min_coord, range);
      }
	};


	// Create the threads
	std::vector<std::thread> workers;
	for (int i = 0; i < num_threads; ++i)
		workers.push_back(std::thread(worker_fn, i));
	for (auto& t : workers)	t.join();
}

void compute_morton_code_openmp(const int input_size, std::vector<Eigen::Vector3f>& inputs,std::vector<Code_t>& morton_keys, float min_coord, const float range, int num_threads ){

	const auto elements_per_thread = math::divide_ceil<int>(input_size, num_threads);
    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        int start = thread_id * elements_per_thread;
        int end = (thread_id == num_threads - 1) ? input_size : (thread_id + 1) * elements_per_thread;

        for (int t = start; t < end; ++t) {
            Eigen::Vector3f input = inputs[t];
            morton_keys[t] = PointToCode(input.x(), input.y(), input.z(), min_coord, range);
        }
    }
}

void compute_morton_code_openmp_pcl(const int input_size, std::vector<pcl::PointXYZ, Eigen::aligned_allocator<pcl::PointXYZ>>& inputs,std::vector<Code_t>& morton_keys, float min_coord, const float range, int num_threads ){

	const auto elements_per_thread = math::divide_ceil<int>(input_size, num_threads);
    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        int start = thread_id * elements_per_thread;
        int end = (thread_id == num_threads - 1) ? input_size : (thread_id + 1) * elements_per_thread;

        for (int t = start; t < end; ++t) {
            pcl::PointXYZ input = inputs[t];
            morton_keys[t] = PointToCode(input.x, input.y, input.z, min_coord, range);
        }
    }
}



void PrefixSumThreaded(const std::vector<int>& edge_counts, std::vector<int>& oc_node_offsets, int n, int num_threads){
    int *suma;
    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
        const int ithread = omp_get_thread_num();
        const int nthreads = omp_get_num_threads();
        #pragma omp single
        {
            suma = new int[nthreads+1];
            suma[0] = 0;
        }
        int sum = 0;
        #pragma omp for schedule(static)
        for (int i=0; i<n; i++) {
            sum += edge_counts[i];
            oc_node_offsets[i+1] = sum;
        }
        suma[ithread+1] = sum;
        #pragma omp barrier
        int offset = 0;
        for(int i=0; i<(ithread+1); i++) {
            offset += suma[i];
        }
        #pragma omp for schedule(static)
        for (int i=0; i<n; i++) {
            oc_node_offsets[i+1] += offset;
        }
    }
    delete[] suma;
}
