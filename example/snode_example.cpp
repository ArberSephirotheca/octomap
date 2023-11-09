#include "runtime/llvm/llvm_runtime_executor.h"
#include "ir/snode.h"
using namespace redwood::lang;
void run_snode(){
    auto compile_config = CompileConfig();
    auto executor = LlvmRuntimeExecutor(compile_config);
    int n = 10;
    executor.materialize_runtime();
    auto *root = new SNode(0, SNodeType::root);
    auto *pointer = &root->pointer(Axis(0), n);
    auto 
}