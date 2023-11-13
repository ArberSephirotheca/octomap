#include "runtime/llvm/llvm_struct_compiler.h"

namespace redwood::lang {
    LlvmStructCompiler::LlvmStructCompiler(Arch arch, const CompileConfig &config, int snode_tree_id)
        : arch_(arch), snode_tree_id_(snode_tree_id) {
    }
    
    void LlvmStructCompiler::run(SNode &root){
        collect_snodes(root);
        auto snodes_rev = snodes;
        std::reverse(snodes_rev.begin(), snodes_rev.end());
    }
}