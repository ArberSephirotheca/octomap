#pragma once
#include "struct/struct.h"
#include "struct/snode.h"
#include "common/arch.h"
#include "program/compile_config.h"

namespace redwood::lang {
class vulkanStructCompiler : public StructCompiler{
    public:
    vulkanStructCompiler(Arch arch, const CompileConfig &config, int snode_tree_id);
    void run(SNode &node) override;

    private:
    Arch arch_;
    int snode_tree_id_;
};
}