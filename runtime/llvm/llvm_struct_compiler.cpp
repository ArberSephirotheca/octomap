#include "runtime/llvm/llvm_struct_compiler.h"

namespace redwood::lang {

    void LlvmStructCompiler::run(SNode &node){
        collect_snodes(root);
        auto snodes_rev = snodes;
        std::reverse(snodes_rev.begin(), snodes_rev.end());
    }
}