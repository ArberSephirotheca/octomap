#pragma once

#include "redwood/codegen/compiled_kernel_data.h"
#include "redwood/program/launch_context_builder.h"

namespace redwood::lang {

class KernelLauncher {
 public:
  using Handle = KernelLaunchHandle;

  virtual void launch_kernel(const CompiledKernelData &compiled_kernel_data,
                             LaunchContextBuilder &ctx) = 0;

  virtual ~KernelLauncher() = default;
};

}  // namespace redwood::lang