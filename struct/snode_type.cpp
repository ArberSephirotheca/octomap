#include "struct/snode_types.h"

#include "common/logging.h"

namespace redwood::lang {

std::string snode_type_name(SNodeType t) {
  switch (t) {
#define PER_SNODE(i) \
  case SNodeType::i: \
    return #i;

#include "inc/snodes.inc.h"

#undef PER_SNODE
    default:
      RW_NOT_IMPLEMENTED;
  }
}

bool is_gc_able(SNodeType t) {
  return (t == SNodeType::pointer || t == SNodeType::dynamic);
}

}  // namespace redwood::lang