#pragma once

#include <string>

namespace redwood::lang {

enum class SNodeType {
#define PER_SNODE(x) x,
#include "inc/snodes.inc.h"
#undef PER_SNODE
};

std::string snode_type_name(SNodeType t);

bool is_gc_able(SNodeType t);

}  // namespace redwood::lang