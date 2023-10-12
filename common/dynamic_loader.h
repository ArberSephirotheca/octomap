#pragma once

#include "common/core.h"

namespace redwood {

class DynamicLoader {
 public:
  static bool check_lib_loaded(const std::string &lib_path);

 private:
  void load_dll(const std::string &dll_path);

  void close_dll();

 public:
  explicit DynamicLoader(const std::string &dll_path);

  void *load_function(const std::string &func_name);

  template <typename T>
  void load_function(const std::string &func_name, T &f) {
    f = (T)load_function(func_name);
  }

  bool loaded() const;

  ~DynamicLoader();

 private:
  void *dll_ = nullptr;
};

}  // namespace redwood