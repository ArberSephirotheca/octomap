#pragma once
#include <cstddef>
constexpr int redwood_max_num_indices = 12;
constexpr int redwood_max_num_args_total = 64;
constexpr int redwood_max_num_args_extra = 32;
constexpr int redwood_max_num_snodes = 1024;
constexpr int kMaxNumSnodeTreesLlvm = 512;
constexpr int redwood_max_gpu_block_dim = 1024;
constexpr std::size_t redwood_global_tmp_buffer_size = 1024 * 1024;
constexpr int redwood_max_num_mem_requests = 1024 * 64;
constexpr std::size_t redwood_page_size = 4096;
constexpr std::size_t redwood_error_message_max_length = 2048;
constexpr std::size_t redwood_error_message_max_num_arguments = 32;
constexpr std::size_t redwood_result_buffer_entries = 32;
constexpr std::size_t redwood_max_num_ret_value = 30;

enum class ParameterType {
  kScalar,
  kNdarray,
  kTexture,
  kRWTexture,
  kTensor,
  kUnknown
};

enum class ExternalArrayLayout { kAOS, kSOA, kNull };

enum class AutodiffMode { kForward, kReverse, kNone, kCheckAutodiffValid };

enum class SNodeGradType { kPrimal, kAdjoint, kDual, kAdjointCheckbit };

enum class BoundaryMode { kUnsafe, kClamp };