#include "runtime/llvm/locked_task.h"
int atomic_max_i32(std::atomic<uint64_t>* dest, int val){
  uint64_t old = dest->load();

  if (val > old) {
    dest->compare_exchange_strong(old, val);
  }
  return old;

}
void mutex_lock_i32(std::atomic<uint64_t>* mutex){
  while(std::atomic_exchange(mutex, 1) == 1);
}
void mutex_unlock_i32(std::atomic<uint64_t>* mutex) {
  std::atomic_exchange(mutex, false);
}

