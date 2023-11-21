#include <limits>
#include <algorithm>
#include <cstring>
#include "struct/snode.h"
#include "common/logging.h"
#include "util/bit.h"
#include "runtime/llvm/runtime.h"
#include "runtime/llvm/node_dynamic.h"
#include "runtime/llvm/node_pointer.h"
/*
#include "redwood/ir/ir.h"
#include "redwood/ir/statements.h"
#include "redwood/program/program.h"
#include "redwood/program/snode_rw_accessors_bank.h"
*/
namespace redwood::lang {

std::atomic<int> SNode::counter{0};

using Ptr = uint8_t *;


Ptr SNode::operator[](int index){
  Ptr ptr = nullptr;
  if (this->type == SNodeType::place){
    // get the parent node, that is where the data is stored
    auto parent = this->parent;
    if (parent->type == SNodeType::pointer){
      auto pointer = static_cast<PointerNode*>(parent);
      ptr =  pointer->Pointer_lookup_element(index);
      // not actiavte yet(or allocating yet, allocate first and then lookup again)
      if (ptr == nullptr){
        pointer->Pointer_activate(index);
        ptr =  pointer->Pointer_lookup_element(index);
    }
    }
  } else if (this->type == SNodeType::root){
    RW_ERROR("SNode type is root, cannot access data.");
  }
  // TODO: add
  else{
    RW_ERROR("SNode type is not place, cannot access data.");
  }
  return ptr;
}


SNode &SNode::insert_children(SNodeType t) {
  RW_ASSERT(t != SNodeType::root);
  std::unique_ptr<SNode> new_ch{nullptr};
   if (t == SNodeType::dynamic){
    new_ch = std::make_unique<DynamicNode>(depth + 1, t);
  }else{
    new_ch = std::make_unique<SNode>(depth + 1, t);
  }
  new_ch->parent = this;
  new_ch->is_path_all_dense = (is_path_all_dense && !new_ch->need_activation());
  for (int i = 0; i < redwood_max_num_indices; i++) {
    new_ch->extractors[i].num_elements_from_root *=
        extractors[i].num_elements_from_root;
  }
  std::memcpy(new_ch->physical_index_position, physical_index_position,
              sizeof(physical_index_position));
  new_ch->num_active_indices = num_active_indices;
  if (type == SNodeType::bit_struct || type == SNodeType::quant_array) {
    new_ch->is_bit_level = true;
  } else {
    new_ch->is_bit_level = is_bit_level;
  }
  ch.push_back(std::move(new_ch));
  return *ch.back();
}

SNode &SNode::create_node(std::vector<Axis> axes,
                          std::vector<int> sizes,
                          SNodeType type,
                          redwood::runtime::LLVMRuntime *runtime
                          /**const DebugInfo &dbg_info*/) {
  if (sizes.size() == 1) {
    sizes = std::vector<int>(axes.size(), sizes[0]);
  }
  if (axes.size() != sizes.size()) {
    RW_ERROR("axes and sizes must have the same size, but got {} and {}.",
             axes.size(), sizes.size());
             /*
    ErrorEmitter(
        redwoodRuntimeError(), &dbg_info,
        fmt::format(
            "axes and sizes must have the same size, but got {} and {}.",
            axes.size(), sizes.size()));
            */
  }

  auto &new_node = insert_children(type);
  for (int i = 0; i < (int)axes.size(); i++) {
    if (sizes[i] <= 0) {
      RW_ERROR("Every dimension of a redwood field should be "
                               "positive, got {} in demension {}.",
                               sizes[i], i);
      
      /*
      ErrorEmitter(redwoodRuntimeError(), &dbg_info,
                   fmt::format("Every dimension of a redwood field should be "
                               "positive, got {} in demension {}.",
                               sizes[i], i));
                               */
    }

    int ind = axes[i].value;
    auto end = new_node.physical_index_position + new_node.num_active_indices;
    bool is_first_division =
        std::find(new_node.physical_index_position, end, ind) == end;
    if (is_first_division) {
      new_node.physical_index_position[new_node.num_active_indices++] = ind;
    } else if (!is_power_of_two(sizes[i])) {
      RW_ERROR(
          "Shape {} is detected on non-first division of axis {}. For "
          "best performance, we recommend that you set it to a power of "
          "two.",
          sizes[i], char('i' + ind));
      /*
      ErrorEmitter(
          redwoodRuntimeWarning(), &dbg_info,
          fmt::format(
              "Shape {} is detected on non-first division of axis {}. For "
              "best performance, we recommend that you set it to a power of "
              "two.",
              sizes[i], char('i' + ind)));
              */
    }
    new_node.extractors[ind].active = true;
    new_node.extractors[ind].num_elements_from_root *= sizes[i];
    new_node.extractors[ind].shape = sizes[i];
  }
  std::sort(new_node.physical_index_position,
            new_node.physical_index_position + new_node.num_active_indices);
  // infer extractors
  int64_t acc_shape = 1;
  for (int i = redwood_max_num_indices - 1; i >= 0; i--) {
    // casting to int32 in extractors.
    new_node.extractors[i].acc_shape = static_cast<int>(acc_shape);
    acc_shape *= new_node.extractors[i].shape;
  }
  if (acc_shape > std::numeric_limits<int>::max()) {
    RW_ERROR("SNode index might be out of int32 boundary but int64 indexing is not "
                 "supported yet. Struct fors might not work either.");
    /*
    ErrorEmitter(
        redwoodIndexWarning(), &dbg_info,
        "SNode index might be out of int32 boundary but int64 indexing is not "
        "supported yet. Struct fors might not work either.");
        */
  }
  new_node.num_cells_per_container = acc_shape;

  if (new_node.type == SNodeType::dynamic) {
    int active_extractor_counder = 0;
    for (int i = 0; i < redwood_max_num_indices; i++) {
      if (new_node.extractors[i].active) {
        active_extractor_counder += 1;
        SNode *p = new_node.parent;
        while (p) {
          if (p->extractors[i].active) {
            RW_ERROR("Dynamic SNode must have a standalone dimensionality.");
            /*
            ErrorEmitter(
                redwoodRuntimeError(), &dbg_info,
                "Dynamic SNode must have a standalone dimensionality.");
                */
          }
          p = p->parent;
        }
      }
    }
    if (active_extractor_counder != 1) {
      RW_ERROR("Dynamic SNode can have only one index extractor.");
      /*
      ErrorEmitter(redwoodRuntimeError(), &dbg_info,
                   "Dynamic SNode can have only one index extractor.");
                   */
    }
  }
  //RW_ASSERT(runtime != nullptr);
  new_node.runtime = runtime;
  return new_node;
}

SNode &SNode::dynamic(const Axis &expr,
                      int n,
                      int chunk_size,
                      redwood::runtime::LLVMRuntime* runtime
                      /*const DebugInfo &dbg_info*/) {
  auto &snode = create_node({expr}, {n}, SNodeType::dynamic/*, dbg_info*/, runtime);
  snode.chunk_size = chunk_size;
  return snode;
}

/*
SNode &SNode::bit_struct(BitStructType *bit_struct_type,
                         const DebugInfo &dbg_info) {
  auto &snode = create_node({}, {}, SNodeType::bit_struct, dbg_info);
  snode.dt = bit_struct_type;
  snode.physical_type = bit_struct_type->get_physical_type();
  return snode;
}
*/


SNode &SNode::quant_array(const std::vector<Axis> &axes,
                          const std::vector<int> &sizes,
                          int bits,
                          redwood::runtime::LLVMRuntime* runtime
                          /*const DebugInfo &dbg_info*/) {
  auto &snode = create_node(axes, sizes, SNodeType::quant_array/*, dbg_info*/, runtime);
  /*
  snode.physical_type =
      TypeFactory::get_instance().get_primitive_int_type(bits, false);
  */
  return snode;
}

bool SNode::is_place() const {
  return type == SNodeType::place;
}

bool SNode::is_scalar() const {
  return is_place() && (num_active_indices == 0);
}

SNode *SNode::get_least_sparse_ancestor() const {
  if (is_path_all_dense) {
    return nullptr;
  }
  auto *result = const_cast<SNode *>(this);
  while (!result->need_activation()) {
    result = result->parent;
    RW_ASSERT(result);
  }
  return result;
}

int SNode::shape_along_axis(int i) const {
  const auto &extractor = extractors[physical_index_position[i]];
  return extractor.num_elements_from_root;
}

/*
int64_t SNode::read_int(const std::vector<int> &i) {
  return snode_rw_accessors_bank_->get(this).read_int(i);
}

uint64_t SNode::read_uint(const std::vector<int> &i) {
  return snode_rw_accessors_bank_->get(this).read_uint(i);
}

float SNode::read_float(const std::vector<int> &i) {
  return snode_rw_accessors_bank_->get(this).read_float(i);
}

void SNode::write_int(const std::vector<int> &i, int64_t val) {
  snode_rw_accessors_bank_->get(this).write_int(i, val);
}

void SNode::write_uint(const std::vector<int> &i, uint64_t val) {
  snode_rw_accessors_bank_->get(this).write_uint(i, val);
}

void SNode::write_float(const std::vector<int> &i, float val) {
  snode_rw_accessors_bank_->get(this).write_float(i, val);
}

*/


SNode::SNode(int depth,
             SNodeType t)
    : depth(depth),
      type(t) {
  id = counter++;
  //node_type_name = get_node_type_name();
  num_active_indices = 0;
  std::memset(physical_index_position, -1, sizeof(physical_index_position));
  parent = nullptr;
  has_ambient = false;
  //dt = PrimitiveType::gen;
  _morton = false;
}

SNode::SNode(const SNode &) {
  RW_NOT_IMPLEMENTED;  // Copying an SNode is forbidden. However we need the
                       // definition here to make pybind11 happy.
}

/*
std::string SNode::get_node_type_name() const {
  return fmt::format("S{}", id);
}
*/

/*
std::string SNode::get_node_type_name_hinted() const {
  std::string suffix;
  if (type == SNodeType::place || type == SNodeType::bit_struct)
    suffix = fmt::format("<{}>", dt->to_string());
  if (is_bit_level)
    suffix += "<bit>";
  return fmt::format("S{}{}{}", id, snode_type_name(type), suffix);
}


void SNode::print() {
  for (int i = 0; i < depth; i++) {
    fmt::print("  ");
  }
  fmt::print("{}", get_node_type_name_hinted());
  fmt::print("\n");
  for (auto &c : ch) {
    c->print();
  }
}
*/

void SNode::set_index_offsets(std::vector<int> index_offsets_) {
  RW_ASSERT(this->index_offsets.empty());
  RW_ASSERT(!index_offsets_.empty());
  RW_ASSERT(type == SNodeType::place);
  RW_ASSERT(index_offsets_.size() == this->num_active_indices);
  this->index_offsets = index_offsets_;
}

// TODO: rename to is_sparse?
bool SNode::need_activation() const {
  return type == SNodeType::pointer || type == SNodeType::hash ||
         type == SNodeType::bitmasked || type == SNodeType::dynamic;
}

/*
void SNode::lazy_grad() {
  make_lazy_place(
      this, snode_to_fields_,
      [this](std::unique_ptr<SNode> &c, std::vector<Expr> &new_grads) {
        if (c->type == SNodeType::place && c->is_primal() && is_real(c->dt) &&
            !c->has_adjoint()) {
          new_grads.push_back(snode_to_fields_->at(c.get())->adjoint);
        }
      });
}

void SNode::lazy_dual() {
  make_lazy_place(
      this, snode_to_fields_,
      [this](std::unique_ptr<SNode> &c, std::vector<Expr> &new_duals) {
        if (c->type == SNodeType::place && c->is_primal() && is_real(c->dt) &&
            !c->has_dual()) {
          new_duals.push_back(snode_to_fields_->at(c.get())->dual);
        }
      });
}

void SNode::allocate_adjoint_checkbit() {
  make_lazy_place(this, snode_to_fields_,
                  [this](std::unique_ptr<SNode> &c,
                         std::vector<Expr> &new_adjoint_checkbits) {
                    if (c->type == SNodeType::place && c->is_primal() &&
                        is_real(c->dt) && c->has_adjoint()) {
                      new_adjoint_checkbits.push_back(
                          snode_to_fields_->at(c.get())->adjoint_checkbit);
                    }
                  });
}
*/

bool SNode::is_primal() const {
  return grad_info && grad_info->is_primal();
}

SNodeGradType SNode::get_snode_grad_type() const {
  RW_ASSERT(grad_info);
  return grad_info->get_snode_grad_type();
}

bool SNode::has_adjoint() const {
  return is_primal() && (grad_info->adjoint_snode() != nullptr);
}

bool SNode::has_adjoint_checkbit() const {
  return is_primal() && (grad_info->adjoint_checkbit_snode() != nullptr);
}

bool SNode::has_dual() const {
  return is_primal() && (grad_info->dual_snode() != nullptr);
}

SNode *SNode::get_adjoint() const {
  RW_ASSERT(has_adjoint());
  return grad_info->adjoint_snode();
}

SNode *SNode::get_adjoint_checkbit() const {
  // RW_ASSERT(has_adjoint());
  return grad_info->adjoint_checkbit_snode();
}

SNode *SNode::get_dual() const {
  RW_ASSERT(has_dual());
  return grad_info->dual_snode();
}

void SNode::set_snode_tree_id(int id) {
  snode_tree_id_ = id;
  for (auto &child : ch) {
    child->set_snode_tree_id(id);
  }
}

int SNode::get_snode_tree_id() const {
  return snode_tree_id_;
}

const SNode *SNode::get_root() const {
  if (!parent) {  // root->parent == nullptr
    return this;
  }
  return parent->get_root();
}

}  // namespace redwood::lang