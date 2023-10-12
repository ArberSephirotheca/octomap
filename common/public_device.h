#pragma once

#include <string>
#include <vector>
#include <assert.h>
#include <memory>

namespace redwood::lang{
class Device;
struct DeviceAllocation;
struct DevicePtr;
using DeviceAllocationId = uint64_t;

enum class PipelineSourceType {
  spirv_binary,
  spirv_src,
  glsl_src,
  hlsl_src,
  dxil_binary,
  llvm_ir_src,
  llvm_ir_binary,
  metal_src,
  metal_ir
};

enum class PipelineStageType {
  compute,
  fragment,
  vertex,
  tesselation_control,
  tesselation_eval,
  geometry,
  raytracing
};

struct PipelineSourceDesc {
  PipelineSourceType type;
  const void *data{nullptr};
  size_t size{0};
  PipelineStageType stage{PipelineStageType::compute};
};

enum class RedwoodResult {
  success = 0,
  error = -1,
  invalid_usage = -2,
  not_supported = -3,
  out_of_memory = -4,
};

#define MAKE_ENUM_FLAGS(name)                  \
  inline name operator|(name a, name b) {      \
    return static_cast<name>(int(a) | int(b)); \
  }                                            \
  inline name operator&(name a, name b) {      \
    return static_cast<name>(int(a) & int(b)); \
  }                                            \
  inline bool operator&&(name a, name b) { return (int(a) & int(b)) != 0; }

enum class AllocUsage : int {
  None = 0,
  Storage = 1,
  Uniform = 2,
  Vertex = 4,
  Index = 8,
  Upload = 16,
};


struct  DeviceAllocation {
  Device *device{nullptr};
  DeviceAllocationId alloc_id{0};
  // TODO: Shall we include size here?

  DevicePtr get_ptr(uint64_t offset = 0) const;

  bool operator==(const DeviceAllocation &other) const {
    return other.device == device && other.alloc_id == alloc_id;
  }

  bool operator!=(const DeviceAllocation &other) const {
    return !(*this == other);
  }
};



struct  DeviceAllocationGuard : public DeviceAllocation {
  explicit DeviceAllocationGuard(DeviceAllocation alloc)
      : DeviceAllocation(alloc) {
  }
  DeviceAllocationGuard(const DeviceAllocationGuard &) = delete;
  ~DeviceAllocationGuard();
};

using DeviceAllocationUnique = std::unique_ptr<DeviceAllocationGuard>;


struct  DevicePtr : public DeviceAllocation {
  uint64_t offset{0};

  bool operator==(const DevicePtr &other) const {
    return other.device == device && other.alloc_id == alloc_id &&
           other.offset == offset;
  }

  bool operator!=(const DevicePtr &other) const {
    return !(*this == other);
  }
};

constexpr DeviceAllocation kDeviceNullAllocation{};
constexpr DevicePtr kDeviceNullPtr{};

class Pipeline {
 public:
  virtual ~Pipeline() {
  }
};

using UPipeline = std::unique_ptr<Pipeline>;


MAKE_ENUM_FLAGS(AllocUsage)

class  PipelineCache {
 public:
  virtual ~PipelineCache() = default;

  /**
   * Get the pointer to the raw data of the cache.
   * - Can return `nullptr` if cache is invalid or empty.
   */
  virtual void *data() noexcept {
    return nullptr;
  }

  /**
   * Get the size of the cache (in bytes).
   */
  virtual size_t size() const noexcept {
    return 0;
  }
};

using UPipelineCache = std::unique_ptr<PipelineCache>;

class CommandList {
 public:
  virtual ~CommandList() {
  }

  /**
   * Bind a pipeline to the command list.
   * Doing so resets all bound resources.
   * @params[in] pipeline The pipeline to be bound
   */
  virtual void bind_pipeline(Pipeline *p) noexcept = 0;

  /**
   * Insert a memory barrier into the command list.
   * The barrier affects a continous region of memory.
   * Changes to memory before the barrier will be visible to accesses after the
   * barrier (API command ordering). i.e. Command later to this barrier will see
   * the changes made by commands before this barrier.
   * This barrier is limited in scope to the Stream that the command list is
   * submitted to. Other Streams or Devices may not observe this barrier.
   * @params[in] ptr The pointer to the start of the region
   * @params[in] size The size of the memory region.
   *                  Size is clamped to the underlying buffer size.
   */
  virtual void buffer_barrier(DevicePtr ptr, size_t size) noexcept = 0;

  /**
   * Insert a memory barrier into the command list.
   * The barrier affects an entire buffer.
   * Behaviour is the same as `buffer_barrier(DevicePtr, size_t)`
   * @params[in] alloc The memory allocation of this barrier
   */
  virtual void buffer_barrier(DeviceAllocation alloc) noexcept = 0;

  /**
   * Insert a memory barrier into the command list.
   * The barrier affects all global memory.
   * Behaviour is the same as `buffer_barrier(DevicePtr, size_t)`
   * @params[in] alloc The memory allocation of this barrier
   */
  virtual void memory_barrier() noexcept = 0;

  /**
   * Insert a buffer copy operation into the command list.
   * @params[in] src The source Device Pointer
   * @params[in] dst The destination Device Pointer
   * @params[in] size The size of the region to be copied.
   *                  The size will be clamped to the minimum between
   *                  `dst.size - dst.offset` and `src.size - src.offset`
   */
  virtual void buffer_copy(DevicePtr dst,
                           DevicePtr src,
                           size_t size) noexcept = 0;

  /**
   * Insert a memory region fill operation into the command list
   * The memory region will be filled with the given (bit precise) value.
   * - (Encouraged behavior) If the `data` is 0, the underlying API might
   *   provide a faster code path.
   * - (Encouraged behavior) If the `size` is -1 (max of size_t) the underlying
   *   API might provide a faster code path.
   * @params[in] ptr The start of the memory region.
   * - ptr.offset will be aligned down to a multiple of 4 bytes.
   * @params[in] size The size of the region.
   * - The size will be clamped to the underlying buffer's size.
   */
  virtual void buffer_fill(DevicePtr ptr,
                           size_t size,
                           uint32_t data) noexcept = 0;

  /**
   * Enqueues a compute operation with {X, Y, Z} amount of workgroups.
   * The block size / workgroup size is pre-determined within the pipeline.
   * - This is only valid if the pipeline has a predetermined block size
   * - This API has a device-dependent variable max values for X, Y, Z
   * - The currently bound pipeline will be dispatched
   * - The enqueued operation starts in CommandList API ordering.
   * - The enqueued operation may end out-of-order, but it respects barriers
   * @params[in] x The number of workgroups in X dimension
   * @params[in] y The number of workgroups in Y dimension
   * @params[in] z The number of workgroups in Y dimension
   * @return The status of this operation
   * - `success` if the operation is successful
   * - `invalid_operation` if the current pipeline has variable block size
   * - `not_supported` if the requested X, Y, or Z is not supported
   */
  virtual RedwoodResult dispatch(uint32_t x,
                             uint32_t y = 1,
                             uint32_t z = 1) noexcept = 0;

  struct ComputeSize {
    uint32_t x{0};
    uint32_t y{0};
    uint32_t z{0};
  };

  /**
   * Enqueues a compute operation with `grid_size` amount of threads.
   * The workgroup size is dynamic and specified through `block_size`
   * - This is only valid if the pipeline has a predetermined block size
   * - This API has a device-dependent variable max values for `grid_size`
   * - This API has a device-dependent supported values for `block_size`
   * - The currently bound pipeline will be dispatched
   * - The enqueued operation starts in CommandList API ordering.
   * - The enqueued operation may end out-of-order, but it respects barriers
   * @params[in] grid_size The number of threads dispatch
   * @params[in] block_size The shape of each block / workgroup / threadsgroup
   * @return The status of this operation
   * - `success` if the operation is successful
   * - `invalid_operation` if the current pipeline has variable block size
   * - `not_supported` if the requested sizes are not supported
   * - `error` if the operation failed due to other reasons
   */
  virtual RedwoodResult dispatch(ComputeSize grid_size,
                             ComputeSize block_size) noexcept {
    return RedwoodResult::not_supported;
  }

  // Profiler support
  virtual void begin_profiler_scope(const std::string &kernel_name) {
  }

  virtual void end_profiler_scope() {
  }
};


class StreamSemaphoreObject {
 public:
  virtual ~StreamSemaphoreObject() {
  }
};

using StreamSemaphore = std::shared_ptr<StreamSemaphoreObject>;

class  Stream {
 public:
  virtual ~Stream() {
  }

  /**
   * Allocates a new CommandList object from the stream.
   * @params[out] out_cmdlist The allocated command list.
   * @return The status of this operation.
   * - `success` If allocation succeeded.
   * - `out_of_memory` If allocation failed due to lack of device or host
   * memory.
   */
  virtual RedwoodResult new_command_list(CommandList **out_cmdlist) noexcept = 0;

  inline std::pair<std::unique_ptr<CommandList>, RedwoodResult>
  new_command_list_unique() {
    CommandList *cmdlist{nullptr};
    RedwoodResult res = this->new_command_list(&cmdlist);
    return std::make_pair(std::unique_ptr<CommandList>(cmdlist), res);
  }

  virtual StreamSemaphore submit(
      CommandList *cmdlist,
      const std::vector<StreamSemaphore> &wait_semaphores = {}) = 0;
  virtual StreamSemaphore submit_synced(
      CommandList *cmdlist,
      const std::vector<StreamSemaphore> &wait_semaphores = {}) = 0;

  virtual void command_sync() = 0;
};

class  Device {

 public:
  virtual ~Device(){};

  struct AllocParams {
    uint64_t size{0};
    bool host_write{false};
    bool host_read{false};
    bool export_sharing{false};
    AllocUsage usage{AllocUsage::Storage};
  };

  virtual RedwoodResult allocate_memory(const AllocParams &params,
                                    DeviceAllocation *out_devalloc) = 0;

  virtual void dealloc_memory(DeviceAllocation handle) = 0;

  virtual uint64_t get_memory_physical_pointer(DeviceAllocation handle) {
    // FIXME: (penguinliong) This method reports the actual device memory
    // address, it's used for bindless (like argument buffer on Metal). If the
    // backend doesn't have access to physical memory address, it should return
    // null and it depends on the backend implementation to use the address in
    // argument binders.
    return 0;
  }

  /**
   * Create a Pipeline Cache, which acclerates backend API's pipeline creation.
   * @params[out] out_cache The created pipeline cache.
   * - If operation failed this will be set to `nullptr`
   * @params[in] initial_size Size of the initial data, can be 0.
   * @params[in] initial_data The initial data, can be nullptr.
   * - This data can be used to load back the cache from previous invocations.
   * - The backend API may ignore this data or deem it incompatible.
   * @return The status of this operation.
   * - `success` if the pipeline cache is created successfully.
   * - `out_of_memory` if operation failed due to lack of device or host memory.
   * - `error` if operation failed due to other errors.
   */
  virtual RedwoodResult create_pipeline_cache(
      PipelineCache **out_cache,
      size_t initial_size = 0,
      const void *initial_data = nullptr) noexcept {
    *out_cache = nullptr;
    return RedwoodResult::not_supported;
  }

  inline std::pair<UPipelineCache, RedwoodResult> create_pipeline_cache_unique(
      size_t initial_size = 0,
      const void *initial_data = nullptr) noexcept {
    PipelineCache *cache{nullptr};
    RedwoodResult res =
        this->create_pipeline_cache(&cache, initial_size, initial_data);
    return std::make_pair(UPipelineCache(cache), res);
  }

  /**
   * Create a Pipeline. A Pipeline is a program that can be dispatched into a
   * stream through a command list.
   * @params[out] out_pipeline The created pipeline.
   * @params[in] src The source description of the pipeline.
   * @params[in] name The name of such pipeline, for debug purposes.
   * @params[in] cache The pipeline cache to use, can be nullptr.
   * @return The status of this operation.
   * - `success` if the pipeline is created successfully.
   * - `out_of_memory` if operation failed due to lack of device or host memory.
   * - `invalid_usage` if the specified source is incompatible or invalid.
   * - `not_supported` if the pipeline uses features the device can't support.
   * - `error` if the operation failed due to other reasons.
   */
  virtual RedwoodResult create_pipeline(
      Pipeline **out_pipeline,
      const PipelineSourceDesc &src,
      std::string name = "Pipeline",
      PipelineCache *cache = nullptr) noexcept = 0;

  inline std::pair<UPipeline, RedwoodResult> create_pipeline_unique(
      const PipelineSourceDesc &src,
      std::string name = "Pipeline",
      PipelineCache *cache = nullptr) noexcept {
    Pipeline *pipeline{nullptr};
    RedwoodResult res = this->create_pipeline(&pipeline, src, name, cache);
    return std::make_pair(UPipeline(pipeline), res);
  }

  inline std::pair<DeviceAllocationUnique, RedwoodResult> allocate_memory_unique(
      const AllocParams &params) {
    DeviceAllocation alloc;
    RedwoodResult res = allocate_memory(params, &alloc);
    if (res != RedwoodResult::success) {
      return std::make_pair(nullptr, res);
    }
    return std::make_pair(std::make_unique<DeviceAllocationGuard>(alloc), res);
  }

  /**
   * Upload data to device allocations immediately.
   * - This is a synchronous operation, function returns when upload is complete
   * - The host data pointers must be valid and large enough for the size of the
   * copy, otherwise this function might segfault
   * - `device_ptr`, `data`, and `sizes` must contain `count` number of valid
   * values
   * @params[in] device_ptr The array to destination device pointers.
   * @params[in] data The array to source host pointers.
   * @params[in] sizes The array to sizes of data/copy.
   * @params[in] count The number of uploads to perform.
   * @return The status of this operation
   * - `success` if the upload is successful.
   * - `out_of_memory` if operation failed due to lack of device or host memory.
   * - `invalid_usage` if the specified source is incompatible or invalid.
   * - `error` if the operation failed due to other reasons.
   */
  virtual RedwoodResult upload_data(DevicePtr *device_ptr,
                                const void **data,
                                size_t *size,
                                int num_alloc = 1) noexcept;

  /**
   * Read data from device allocations back to host immediately.
   * - This is a synchronous operation, function returns when readback is
   * complete
   * - The host data pointers must be valid and large enough for the size of the
   * copy, otherwise this function might segfault
   * - `device_ptr`, `data`, and `sizes` must contain `count` number of valid
   * values
   * @params[in] device_ptr The array to source device pointers.
   * @params[in] data The array to destination host pointers.
   * @params[in] sizes The array to sizes of data/copy.
   * @params[in] count The number of readbacks to perform.
   * @params[in] wait_sema The semaphores to wait for before the copy is
   * initiated.
   * @return The status of this operation
   * - `success` if the upload is successful.
   * - `out_of_memory` if operation failed due to lack of device or host memory.
   * - `invalid_usage` if the specified source is incompatible or invalid.
   * - `error` if the operation failed due to other reasons.
   */
  virtual RedwoodResult readback_data(
      DevicePtr *device_ptr,
      void **data,
      size_t *size,
      int num_alloc = 1,
      const std::vector<StreamSemaphore> &wait_sema = {}) noexcept;

  // Each thraed will acquire its own stream
  virtual Stream *get_compute_stream() = 0;

  // Wait for all tasks to complete (task from all streams)
  virtual void wait_idle() = 0;

  /**
   * Map a range within a DeviceAllocation memory into host address space.
   *
   * @param[in] ptr The Device Pointer to map.
   * @param[in] size The size of the mapped region.
   * @param[out] mapped_ptr Outputs the pointer to the mapped region.
   * @return The result status.
   *         `success` when the mapping is successful.
   *         `invalid_usage` when the memory is not host visible.
   *         `invalid_usage` when trying to map the memory multiple times.
   *         `invalid_usage` when `ptr.offset + size` is out-of-bounds.
   *         `error` when the mapping failed for other reasons.
   */
  virtual RedwoodResult map_range(DevicePtr ptr,
                              uint64_t size,
                              void **mapped_ptr) = 0;

  /**
   * Map an entire DeviceAllocation into host address space.
   * @param[in] ptr The Device Pointer to map.
   * @param[in] size The size of the mapped region.
   * @param[out] mapped_ptr Outputs the pointer to the mapped region.
   * @return The result status.
   *         `success` when the mapping is successful.
   *         `invalid_usage` when the memory is not host visible.
   *         `invalid_usage` when trying to map the memory multiple times.
   *         `invalid_usage` when `ptr.offset + size` is out-of-bounds.
   *         `error` when the mapping failed for other reasons.
   */
  virtual RedwoodResult map(DeviceAllocation alloc, void **mapped_ptr) = 0;

  /**
   * Unmap a previously mapped DevicePtr or DeviceAllocation.
   * @param[in] ptr The DevicePtr to unmap.
   */
  virtual void unmap(DevicePtr ptr) = 0;

  /**
   * Unmap a previously mapped DevicePtr or DeviceAllocation.
   * @param[in] alloc The DeviceAllocation to unmap
   */
  virtual void unmap(DeviceAllocation alloc) = 0;

  // Directly share memory in the form of alias
  static DeviceAllocation share_to(DeviceAllocation *alloc, Device *target);

  // Strictly intra device copy (synced)
  virtual void memcpy_internal(DevicePtr dst, DevicePtr src, uint64_t size) = 0;

  // Copy memory inter or intra devices (synced)
  enum class MemcpyCapability { Direct, RequiresStagingBuffer, RequiresHost };

  static MemcpyCapability check_memcpy_capability(DevicePtr dst,
                                                  DevicePtr src,
                                                  uint64_t size);

  static void memcpy_direct(DevicePtr dst, DevicePtr src, uint64_t size);

  static void memcpy_via_staging(DevicePtr dst,
                                 DevicePtr staging,
                                 DevicePtr src,
                                 uint64_t size);

  static void memcpy_via_host(DevicePtr dst,
                              void *host_buffer,
                              DevicePtr src,
                              uint64_t size);


  // Profiler support
  virtual void profiler_sync() {
  }

  virtual size_t profiler_get_sampler_count() {
    return 0;
  }

  virtual std::vector<std::pair<std::string, double>>
  profiler_flush_sampled_time() {
    return std::vector<std::pair<std::string, double>>();
  }
};
}