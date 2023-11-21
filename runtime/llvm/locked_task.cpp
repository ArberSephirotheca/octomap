#include "runtime/llvm/locked_task.h"
int atomic_max_i32(Ptr dest, int val){
  std::atomic<uint64_t>* d = (std::atomic<uint64_t>*)dest;
  uint64_t old = d->load();

  if (val > old) {
    d->compare_exchange_strong(old, val);
  }
  return old;

}
void mutex_lock_i32(Ptr mutex){
  std::atomic<uint64_t>* m = (std::atomic<uint64_t>*)mutex;
  while(std::atomic_exchange(m, 1) == 1);
}
void mutex_unlock_i32(Ptr mutex) {
  std::atomic<uint64_t>* m = (std::atomic<uint64_t>*)mutex;
  std::atomic_exchange(m, false);
}

