#pragma once

#include "common/core.h"

#include <atomic>
#include <condition_variable>
#include <functional>
#include <thread>

namespace redwood {

using RangeForTaskFunc = void(void *, int thread_id, int i);
using ParallelFor = void(int n, int num_threads, void *, RangeForTaskFunc func);

class ThreadPool {
 public:
  std::vector<std::thread> threads;
  std::condition_variable slave_cv;
  std::condition_variable master_cv;
  std::mutex mutex;
  std::atomic<int> task_head;
  int task_tail;
  int running_threads;
  int max_num_threads;
  int desired_num_threads;
  uint64_t timestamp;
  uint64_t last_finished;
  bool started;
  bool exiting;
  RangeForTaskFunc *func;
  void *range_for_task_context;  // Note: this is a pointer to a
                                 // range_task_helper_context defined in the
                                 // LLVM runtime, which is different from
                                 // redwood::lang::Context.
  int thread_counter;

  explicit ThreadPool(int max_num_threads);

  void run(int splits,
           int desired_num_threads,
           void *range_for_task_context,
           RangeForTaskFunc *func);

  static void static_run(ThreadPool *pool,
                         int splits,
                         int desired_num_threads,
                         void *range_for_task_context,
                         RangeForTaskFunc *func) {
    return pool->run(splits, desired_num_threads, range_for_task_context, func);
  }

  void target();

  ~ThreadPool();
};

}  // namespace redwood