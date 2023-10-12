#pragma once

#include "common/core.h"

#include <string>
#include <cstdlib>

namespace redwood::lang {

static inline int get_environ_config(const std::string &name,
                                     int default_value = 0) {
  char *res = std::getenv(name.c_str());
  if (res == nullptr)
    return default_value;
  return std::stoi(res);
}

}  // namespace redwood::lang