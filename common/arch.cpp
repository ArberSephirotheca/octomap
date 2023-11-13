#include <array>
#include <cstdio>
#include <iostream>
#include "common/arch.h"
#include "common/logging.h"
//#include "redwood/rhi/impl_support.h"

namespace redwood {
#define RW_DEBUG_SNPRINTF std::snprintf
#if defined(SPDLOG_H) && defined(RW_WARN)
#define RW_LOG_ERROR(msg) RW_WARN("RHI Error : {}", msg)
#else
#define RW_LOG_ERROR(msg) std::cerr << "RW Error: " << msg << std::endl;
#endif

std::string arch_name(Arch arch) {
  switch (arch) {
#define PER_ARCH(x) \
  case Arch::x:     \
    return #x;      \
    break;
#include "inc/archs.inc.h"

#undef PER_ARCH
    default:
      RW_NOT_IMPLEMENTED
  }
}

Arch arch_from_name(const std::string &arch_name) {
#define PER_ARCH(x)           \
  else if (arch_name == #x) { \
    return Arch::x;           \
  }

  if (false) {
  }
#include "inc/archs.inc.h"

  else {
    std::array<char, 256> buf;
    RW_DEBUG_SNPRINTF(buf.data(), buf.size(), "Unknown architecture name: %s",
                       arch_name.c_str());
    RW_LOG_ERROR(buf.data());
    RW_NOT_IMPLEMENTED
  }

#undef PER_ARCH
}

// Assuming a processor is either a CPU or a GPU. DSP/TPUs not considered.
bool arch_is_cpu(Arch arch) {
  if (arch == Arch::x64 || arch == Arch::arm64 || arch == Arch::js) {
    return true;
  } else {
    return false;
  }
}

bool arch_is_cuda(Arch arch) {
  return arch == Arch::cuda;
}

bool arch_uses_llvm(Arch arch) {
  return (arch == Arch::x64 || arch == Arch::arm64 || arch == Arch::cuda ||
          arch == Arch::dx12 || arch == Arch::amdgpu);
}

bool arch_is_gpu(Arch arch) {
  return !arch_is_cpu(arch);
}

bool arch_uses_spirv(Arch arch) {
  return arch == Arch::opengl || arch == Arch::gles || arch == Arch::vulkan ||
         arch == Arch::dx11 || arch == Arch::metal;
}

Arch host_arch() {
#if defined(RW_ARCH_x64)
  return Arch::x64;
#endif
#if defined(RW_ARCH_ARM)
  return Arch::arm64;
#endif
  RW_NOT_IMPLEMENTED
}

bool arch_use_host_memory(Arch arch) {
  return arch_is_cpu(arch);
}

int default_simd_width(Arch arch) {
  if (arch == Arch::x64) {
    return 8;
  } else if (arch == Arch::cuda) {
    return 32;
  } else if (arch == Arch::arm64) {
    return 4;
  } else {
    RW_NOT_IMPLEMENTED;
    return -1;
  }
}

}  // namespace redwood