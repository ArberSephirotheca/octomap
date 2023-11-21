#pragma once

#include <functional>
#include <cstring>
#include <spdlog/common.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/ostr.h>
namespace spdlog {
class logger;
}

#ifdef _WIN64
#define __FILENAME__ \
  (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#else
#define __FILENAME__ \
  (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

#define SPD_AUGMENTED_LOG(X, ...)                                        \
  redwood::Logger::get_instance().X(                                      \
      fmt::format("[{}:{}@{}] ", __FILENAME__, __FUNCTION__, __LINE__) + \
      fmt::format(__VA_ARGS__))

#if defined(RW_PLATFORM_WINDOWS)
#define RW_UNREACHABLE __assume(0);
#else
#define RW_UNREACHABLE __builtin_unreachable();
#endif

#define RW_TRACE(...) SPD_AUGMENTED_LOG(trace, __VA_ARGS__)
#define RW_DEBUG(...) SPD_AUGMENTED_LOG(debug, __VA_ARGS__)
#define RW_INFO(...) SPD_AUGMENTED_LOG(info, __VA_ARGS__)
#define RW_WARN(...) SPD_AUGMENTED_LOG(warn, __VA_ARGS__)
#define RW_ERROR(...)                      \
  {                                        \
    SPD_AUGMENTED_LOG(error, __VA_ARGS__); \
    RW_UNREACHABLE;                        \
  }
#define RW_CRITICAL(...)                      \
  {                                           \
    SPD_AUGMENTED_LOG(critical, __VA_ARGS__); \
    RW_UNREACHABLE;                           \
  }

#define RW_TRACE_IF(condition, ...) \
  if (condition) {                  \
    RW_TRACE(__VA_ARGS__);          \
  }
#define RW_TRACE_UNLESS(condition, ...) \
  if (!(condition)) {                   \
    RW_TRACE(__VA_ARGS__);              \
  }
#define RW_DEBUG_IF(condition, ...) \
  if (condition) {                  \
    RW_DEBUG(__VA_ARGS__);          \
  }
#define RW_DEBUG_UNLESS(condition, ...) \
  if (!(condition)) {                   \
    RW_DEBUG(__VA_ARGS__);              \
  }
#define RW_INFO_IF(condition, ...) \
  if (condition) {                 \
    RW_INFO(__VA_ARGS__);          \
  }
#define RW_INFO_UNLESS(condition, ...) \
  if (!(condition)) {                  \
    RW_INFO(__VA_ARGS__);              \
  }
#define RW_WARN_IF(condition, ...) \
  if (condition) {                 \
    RW_WARN(__VA_ARGS__);          \
  }
#define RW_WARN_UNLESS(condition, ...) \
  if (!(condition)) {                  \
    RW_WARN(__VA_ARGS__);              \
  }
#define RW_ERROR_IF(condition, ...) \
  if (condition) {                  \
    RW_ERROR(__VA_ARGS__);          \
  }
#define RW_ERROR_UNLESS(condition, ...) \
  if (!(condition)) {                   \
    RW_ERROR(__VA_ARGS__);              \
  }
#define RW_CRITICAL_IF(condition, ...) \
  if (condition) {                     \
    RW_CRITICAL(__VA_ARGS__);          \
  }
#define RW_CRITICAL_UNLESS(condition, ...) \
  if (!(condition)) {                      \
    RW_CRITICAL(__VA_ARGS__);              \
  }

#define RW_ASSERT(x) RW_ASSERT_INFO((x), "Assertion failure: " #x)
#define RW_ASSERT_INFO(x, ...)             \
  {                                        \
    bool ___ret___ = static_cast<bool>(x); \
    if (!___ret___) {                      \
      RW_ERROR(__VA_ARGS__);               \
    }                                      \
  }
#define RW_NOT_IMPLEMENTED RW_ERROR("Not supported.");

#define RW_STOP RW_ERROR("Stopping here")
#define RW_TAG RW_INFO("Tagging here")

#define RW_LOG_SET_PATTERN(x) spdlog::set_pattern(x);

#define RW_FLUSH_LOGGER \
  { redwood::Logger::get_instance().flush(); };

#define RW_P(x) \
  { RW_INFO("{}", redwood::TextSerializer::serialize(#x, (x))); }

namespace redwood {

class Logger {
 private:
  std::shared_ptr<spdlog::logger> console_;
  int level_;
  std::function<void()> print_stacktrace_fn_;

  Logger();

 public:
  void trace(const std::string &s);
  void debug(const std::string &s);
  void info(const std::string &s);
  void warn(const std::string &s);
  void error(const std::string &s, bool raise_exception = true);
  void critical(const std::string &s);
  void flush();
  void set_level(const std::string &level);
  bool is_level_effective(const std::string &level_name);
  int get_level();
  static int level_enum_from_string(const std::string &level);
  void set_level_default();

  // This is mostly to decouple the implementation.
  void set_print_stacktrace_func(std::function<void()> print_fn);

  static Logger &get_instance();
};

}  // namespace redwood