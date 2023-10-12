#include "common/dynamic_loader.h"
#include <dlfcn.h>
namespace redwood {

bool DynamicLoader::check_lib_loaded(const std::string &lib_path) {
  bool loaded = false;
  loaded = (dlopen(lib_path.c_str(), RTLD_NOLOAD) != nullptr);
  return loaded;
}

DynamicLoader::DynamicLoader(const std::string &dll_path) {
  load_dll(dll_path);
}

void DynamicLoader::load_dll(const std::string &dll_path) {
  dll_ = dlopen(dll_path.c_str(), RTLD_LAZY);
}

void *DynamicLoader::load_function(const std::string &func_name) {
  RW_ASSERT_INFO(loaded(), "DLL not opened");
  auto func = dlsym(dll_, func_name.c_str());
  const char *dlsym_error = dlerror();
  RW_ERROR_IF(dlsym_error, "Cannot load function: {}", dlsym_error);
  RW_ERROR_IF(!func, "Function {} not found", func_name);
  return func;
}

void DynamicLoader::close_dll() {
  RW_ASSERT_INFO(loaded(), "DLL not opened");
  dlclose(dll_);
  dll_ = nullptr;
}

DynamicLoader::~DynamicLoader() {
  if (loaded())
    close_dll();
}

bool DynamicLoader::loaded() const {
  return dll_ != nullptr;
}

}  // namespace redwood