#include "common/device.h"

#include "cpu/cpu_device.h"
//#include "cuda/cuda_device.h"

namespace redwood::lang {

DeviceAllocationGuard::~DeviceAllocationGuard() {
  device->dealloc_memory(*this);
}


DevicePtr DeviceAllocation::get_ptr(uint64_t offset) const {
    DevicePtr devicePtr;
    devicePtr.device = this->device;
    devicePtr.alloc_id = this->alloc_id;
    devicePtr.offset = offset;
    return devicePtr;
}

Device::MemcpyCapability Device::check_memcpy_capability(DevicePtr dst,
                                                         DevicePtr src,
                                                         uint64_t size) {
  if (dst.device == src.device) {
    return Device::MemcpyCapability::Direct;
  }

  return Device::MemcpyCapability::RequiresHost;
}

void Device::memcpy_direct(DevicePtr dst, DevicePtr src, uint64_t size) {
  // Intra-device copy
  if (dst.device == src.device) {
    dst.device->memcpy_internal(dst, src, size);
    return;
  }
  RW_NOT_IMPLEMENTED;
}

void Device::memcpy_via_staging(DevicePtr dst,
                                DevicePtr staging,
                                DevicePtr src,
                                uint64_t size) {
  RW_NOT_IMPLEMENTED;
}

void Device::memcpy_via_host(DevicePtr dst,
                             void *host_buffer,
                             DevicePtr src,
                             uint64_t size) {
  RW_NOT_IMPLEMENTED;
}

RedwoodResult Device::upload_data(DevicePtr *device_ptr,
                              const void **data,
                              size_t *size,
                              int num_alloc) noexcept {
  if (!device_ptr || !data || !size) {
    return RedwoodResult::invalid_usage;
  }

  std::vector<DeviceAllocationUnique> stagings;
  for (int i = 0; i < num_alloc; i++) {
    if (device_ptr[i].device != this || !data[i]) {
      return RedwoodResult::invalid_usage;
    }
    auto [staging, res] = this->allocate_memory_unique(
        {size[i], /*host_write=*/true, /*host_read=*/false,
         /*export_sharing=*/false, AllocUsage::Upload});
    if (res != RedwoodResult::success) {
      return res;
    }

    void *mapped{nullptr};
    res = this->map(*staging, &mapped);
    if (res != RedwoodResult::success) {
      return res;
    }
    memcpy(mapped, data[i], size[i]);
    this->unmap(*staging);

    stagings.push_back(std::move(staging));
  }

  Stream *s = this->get_compute_stream();
  auto [cmdlist, res] = s->new_command_list_unique();
  if (res != RedwoodResult::success) {
    return res;
  }
  for (int i = 0; i < num_alloc; i++) {
    cmdlist->buffer_copy(device_ptr[i], stagings[i]->get_ptr(0), size[i]);
  }
  s->submit_synced(cmdlist.get());

  return RedwoodResult::success;
}

RedwoodResult Device::readback_data(
    DevicePtr *device_ptr,
    void **data,
    size_t *size,
    int num_alloc,
    const std::vector<StreamSemaphore> &wait_sema) noexcept {
  if (!device_ptr || !data || !size) {
    return RedwoodResult::invalid_usage;
  }

  Stream *s = this->get_compute_stream();
  auto [cmdlist, res] = s->new_command_list_unique();
  if (res != RedwoodResult::success) {
    return res;
  }

  std::vector<DeviceAllocationUnique> stagings;
  for (int i = 0; i < num_alloc; i++) {
    if (device_ptr[i].device != this || !data[i]) {
      return RedwoodResult::invalid_usage;
    }
    auto [staging, res] = this->allocate_memory_unique(
        {size[i], /*host_write=*/false, /*host_read=*/true,
         /*export_sharing=*/false, AllocUsage::None});
    if (res != RedwoodResult::success) {
      return res;
    }

    cmdlist->buffer_copy(staging->get_ptr(0), device_ptr[i], size[i]);
    stagings.push_back(std::move(staging));
  }
  s->submit_synced(cmdlist.get(), wait_sema);

  for (int i = 0; i < num_alloc; i++) {
    void *mapped{nullptr};
    RedwoodResult res = this->map(*stagings[i], &mapped);
    if (res != RedwoodResult::success) {
      return res;
    }
    memcpy(data[i], mapped, size[i]);
    this->unmap(*stagings[i]);
  }

  return RedwoodResult::success;
}

}  // namespace redwood::lang