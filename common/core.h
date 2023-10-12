#pragma once
#include <iostream>
#include <type_traits>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <string>
#include <functional>
#include "common/platform_macros.h"
#include "common/logging.h"
namespace redwood
{
    template <typename T>
    void trash(T &&t)
    {
        static_assert(!std::is_same<T, void>::value, "");
    }
    class PID
    {
    public:
        static int get_pid();
        static int get_parent_pid();
    };
}