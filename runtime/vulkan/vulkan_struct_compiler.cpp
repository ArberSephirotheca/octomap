#include "runtime/vulkan/vulkan_struct_compiler.h"

namespace redwood::lang {
    vulkanStructCompiler::vulkanStructCompiler(Arch arch, const CompileConfig &config, int snode_tree_id)
        : arch_(arch), snode_tree_id_(snode_tree_id) {
    }
    
    void vulkanStructCompiler::run(SNode &root){
        collect_snodes(root);
        auto snodes_rev = snodes;
        std::reverse(snodes_rev.begin(), snodes_rev.end());
    }
}