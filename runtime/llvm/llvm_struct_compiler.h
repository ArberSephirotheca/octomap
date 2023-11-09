#pragma once
#include "struct/struct.h"
#include "ir/snode.h"

namespace redwood::lang {
class LlvmStructCompiler : public StructCompiler{
    public:
    void run(SNode &node) override;
};
}