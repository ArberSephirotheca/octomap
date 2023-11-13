// Codegen for the hierarchical data structure
#pragma once

#include "struct/snode.h"

namespace redwood::lang {

class StructCompiler {
 public:
  std::vector<SNode *> stack;
  std::vector<SNode *> snodes;
  std::size_t root_size{0};

  virtual ~StructCompiler() = default;

  void collect_snodes(SNode &snode);

  virtual void run(SNode &node) = 0;
};

}  // namespace redwood::lang