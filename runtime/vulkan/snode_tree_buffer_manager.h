#pragma once
#include "inc/constants.h"
#include "struct/snode_tree.h"
#include "common/public_device.h"
#define RW_RUNTIME_HOST

#include <set>
#include <map>

using Ptr = uint8_t *;

namespace redwood::lang {

class JITModule;
class LlvmRuntimeExecutor;

class SNodeTreeBufferManager {
 public:
  explicit SNodeTreeBufferManager(LlvmRuntimeExecutor *runtime_exec);

  Ptr allocate(std::size_t size,
               const int snode_tree_id,
               uint64_t *result_buffer);

  void destroy(SNodeTree *snode_tree);

 private:
  LlvmRuntimeExecutor *runtime_exec_;
  std::map<int, DeviceAllocation> snode_tree_id_to_device_alloc_;
};

}  // namespace redwood::lang