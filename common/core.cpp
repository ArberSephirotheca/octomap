#include "common/core.h"
#include <unistd.h> 

namespace redwood {
int PID::get_pid() {
  return (int)getpid();
}
}