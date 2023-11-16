#include "runtime/llvm/locked_task.h"
int atomic_max_i32(std::atomic<int>& dest, int val){
  int old = dest.load();

  if (val > old) {
    dest.compare_exchange_strong(old, val);
  }
  return old;

}
void mutex_lock_i32(std::atomic<bool>& mutex){
  while(std::atomic_exchange(&mutex, true) == true);
}
void mutex_unlock_i32(std::atomic<bool>& mutex) {
  std::atomic_exchange(&mutex, false);
}

