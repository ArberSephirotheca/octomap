#pragma once

#include <memory>

//#include "vulkan/IR/Module.h"
#include "common/core.h"
#include "struct/snode_types.h"
/*
#include "redwood/common/serialization.h"
#include "redwood/program/kernel.h"
#include "redwood/util/offline_cache.h"
#include "redwood/codegen/vulkan/vulkan_compiled_data.h"
#include "redwood/codegen/vulkan/compiled_kernel_data.h"
*/

namespace redwood::lang {

// NOTE: The vulkanOfflineCache, vulkanOfflineCacheFileReader and
// vulkanOfflineCacheFileWriter are only used by vulkan AOT now.
// TODO(PGZXB): Rename these structs/classes.

struct vulkanOfflineCache {
  using Version = uint16_t[3];  // {MAJOR, MINOR, PATCH}

  enum Format {
    LL = 0x01,
    BC = 0x10,
  };
    /*
  struct KernelCacheData {
    std::string kernel_key;
    std::vector<std::pair<std::vector<int>, Callable::Parameter>> args;
    std::vector<Callable::Ret> rets;
    vulkanCompiledKernel compiled_data;

    const StructType *ret_type = nullptr;
    size_t ret_size{0};

    const StructType *args_type = nullptr;
    size_t args_size{0};

    // For cache cleaning
    std::size_t size{0};          // byte
    std::time_t created_at{0};    // millsec
    std::time_t last_used_at{0};  // millsec

    KernelCacheData() = default;
    KernelCacheData(KernelCacheData &&) = default;
    KernelCacheData &operator=(KernelCacheData &&) = default;
    ~KernelCacheData() = default;

    KernelCacheData clone() const;
    vulkan::CompiledKernelData::InternalData convert_to_vulkan_ckd_data() const;

    TI_IO_DEF(kernel_key,
              args,
              compiled_data,
              size,
              created_at,
              last_used_at,
              rets,
              ret_type,
              ret_size,
              args_type,
              args_size);
  };
  */

  struct FieldCacheData {
    struct SNodeCacheData {
      int id{0};
      SNodeType type = SNodeType::undefined;
      size_t cell_size_bytes{0};
      size_t chunk_size{0};

    };

    int tree_id{0};
    int root_id{0};
    size_t root_size{0};
    std::vector<SNodeCacheData> snode_metas;


    // TODO(zhanlue): refactor vulkan::Modules
    //
    // struct_module will eventually get cloned into each kernel_module,
    // so there's no need to serialize it here.
    //
    // We have three different types of vulkan::Module
    // 1. runtime_module: contains runtime functions.
    // 2. struct_module: contains compiled SNodeTree in vulkan::Type.
    // 3. kernel_modules: contains compiled kernel codes.
    //
    // The way those modules work rely on a recursive clone mechanism:
    //   runtime_module = load("runtime.bc")
    //   struct_module = clone(runtime_module) + compiled-SNodeTree
    //   kernel_module = clone(struct_module) + compiled-Kernel
    //
    // As a result, every kernel_module contains a copy of struct_module +
    // runtime_module.
    //
    // This recursive clone mechanism is super fragile,
    // which potentially causes inconsistency between modules if not handled
    // properly.
    //
    // Let's turn to use vulkan::link to bind the modules,
    // and make runtime_module, struct_module, kernel_module independent of each
    // other
  };
    /*

  using KernelMetadata = KernelCacheData;  // Required by CacheCleaner

  Version version{};
  std::size_t size{0};  // byte

  // TODO(zhanlue): we need a better identifier for each FieldCacheData
  // (SNodeTree) Given that snode_tree_id is not continuous, it is ridiculous to
  // ask the users to remember each of the snode_tree_ids
  // ** Find a way to name each SNodeTree **
  */
  std::unordered_map<int, FieldCacheData> fields;  // key = snode_tree_id
  /*
  std::unordered_map<std::string, KernelCacheData>
      kernels;  // key = kernel_name

  // NOTE: The "version" must be the first field to be serialized
  TI_IO_DEF(version, size, fields, kernels);
  */
};

/*

class vulkanOfflineCacheFileReader {
 public:
  bool get_kernel_cache(vulkanOfflineCache::KernelCacheData &res,
                        const std::string &key,
                        vulkan::vulkanContext &vulkan_ctx);

  bool get_field_cache(vulkanOfflineCache::FieldCacheData &res,
                       int snode_tree_id);

  size_t get_num_snode_trees();

  static std::unique_ptr<vulkanOfflineCacheFileReader> make(
      const std::string &path,
      vulkanOfflineCache::Format format = vulkanOfflineCache::Format::LL);

  static bool load_meta_data(vulkanOfflineCache &data,
                             const std::string &cache_file_path,
                             bool with_lock = true);

 private:
  vulkanOfflineCacheFileReader(const std::string &path,
                             vulkanOfflineCache &&data,
                             vulkanOfflineCache::Format format);

  std::unique_ptr<vulkan::Module> load_module(const std::string &path_prefix,
                                            const std::string &key,
                                            vulkan::vulkanContext &vulkan_ctx) const;

  std::string path_;
  vulkanOfflineCache data_;
  vulkanOfflineCache::Format format_;
};

class vulkanOfflineCacheFileWriter {
 public:
  using CleanCachePolicy = offline_cache::CleanCachePolicy;

  void set_data(vulkanOfflineCache &&data) {
    this->mangled_ = false;
    this->data_ = std::move(data);
  }

  void set_data(std::unique_ptr<vulkanOfflineCache> &&data_ptr) {
    set_data(std::move(*data_ptr.get()));
  }

  void add_kernel_cache(const std::string &key,
                        vulkanOfflineCache::KernelCacheData &&kernel_cache) {
    data_.kernels[key] = std::move(kernel_cache);
  }

  void dump(const std::string &path,
            vulkanOfflineCache::Format format = vulkanOfflineCache::Format::LL,
            bool merge_with_old = false);

  void set_no_mangle() {
    mangled_ = true;
  }

  static void clean_cache(const std::string &path,
                          CleanCachePolicy policy,
                          int max_bytes,
                          double cleaning_factor);

 private:
  void merge_with(vulkanOfflineCache &&data);

  void mangle_offloaded_task_name(const std::string &kernel_key,
                                  vulkanCompiledKernel &compiled_data);

  vulkanOfflineCache data_;
  bool mangled_{false};
};
*/

}  // namespace redwood::lang
