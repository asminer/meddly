
// $Id$

/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2009, Iowa State University Research Foundation, Inc.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.
*/



#include "../include/meddly_expert.h"
#include "../src/operation_ext.h"
#include "../src/compute_cache.h"

#include "config.h"
#include "revision.h"

// ------------------ compute_manager --------------------------
compute_manager::compute_manager() {}


compute_manager::~compute_manager() {}


const char*
compute_manager::getErrorCodeName(compute_manager::error e)
{
  switch (e) {
    case SUCCESS:
        return "Operation returned successfully";
    case NOT_IMPLEMENTED:
        return "Operation not implemented";
    case INSUFFICIENT_MEMORY:
        return "Operation failed -- lack of memory";
    case UNKNOWN_OPERATION:
        return "Operation failed -- bad operation handle";
    case FOREST_MISMATCH:
        return "Operation requires same forest, it was not";
    case TYPE_MISMATCH:
        return
          "Operation not supported for given forest, operand or result type";
    case WRONG_NUMBER:
        return "Operation failed -- incorrect number of operands";
    case OVERFLOW:
        return "Operation failed -- calculation will not fit data type";
    default:
        return "Unknown error code";
  }
}


compute_manager* MEDDLY_getComputeManager()
{
  static expert_compute_manager* ecm = new expert_compute_manager();
  return ecm;
}


const char* MEDDLY_getLibraryInfo(int what)
{
  static char* title = 0;
  switch (what) {
    case 0:
      if (!title) {
        title = new char[80];
        snprintf(title, 80, 
          "%s version %s.%d (32-bit and 64-bit compatible)", 
          PACKAGE_NAME, VERSION, REVISION_NUMBER
        );
      }
      return title;

    case 1:
      return "Copyright (C) 2009, Andrew Miner and Junaid Babar.";

    case 2:
      return "Released under the GNU Lesser General Public License, version 3";
 
    case 3:
      return "http://sourceforge.net/projects/meddly/";

    case 4:
      return "Data Structures and operations available:\n\
(1) MDDs: Union, Intersection, Difference.\n\
(2) Matrix Diagrams (MXDs): Union, Intersection, Difference.\n\
(3) Multi-Terminal MDDs (MTMDDs) with integer or real terminals:\n\
    Arithmetic: Plus, Minus, Multiply, Divide, Min, Max.\n\
    Logical: <, <=, >, >=, ==, !=.\n\
    Conversion to and from MDDs.\n\
(4) Multi-Terminal MXDs (MTMXDs) with integer or real terminals:\n\
    Arithmetic: Plus, Minus, Multiply, Divide, Min, Max.\n\
    Logical: <, <=, >, >=, ==, !=.\n\
    Conversion to and from MXDs.\n\
";
  }
  return 0;
}



// ------------------------ expert_compute_manager --------------------

expert_compute_manager::expert_compute_manager()
{
  // initialize compute cache
  cc = new compute_cache();
  // initialize builtinOpEntries and customOpEntries
  builtinOpEntries = new std::map<builtin_op_key, op_info>();
  customOpEntries = new std::map<custom_op_key, op_info>();
}


expert_compute_manager::~expert_compute_manager()
{
  delete customOpEntries;
  delete builtinOpEntries;
  delete cc;
}


const char* expert_compute_manager::getOperationName(
    compute_manager::op_code op) const
{
  switch(op) {
    case COPY: return "Copy Edge";
    case UNION: return "Union";
    case INTERSECTION: return "Intersection";
    case DIFFERENCE: return "Difference";
    case COMPLEMENT: return "Complement";
    case PRE_IMAGE: return "Pre-Image";
    case POST_IMAGE: return "Post-Image";
    case REACHABLE_STATES_DFS:
                     return "Reachable States via Depth-First Search";
    case REACHABLE_STATES_BFS:
                     return "Reachable States via Breadth-First Search";
    default: return "Unknown operation";
  }
}


compute_manager::error expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, dd_edge &b)
{
  // type check
  compute_manager::error err = owner->op->typeCheck(owner);
  return err != compute_manager::SUCCESS?
    err:
    owner->op->compute(owner, a, b);
}


compute_manager::error expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  // type check
  compute_manager::error err = owner->op->typeCheck(owner);
  return err != compute_manager::SUCCESS?
    err:
    owner->op->compute(owner, a, b, c);
}


compute_manager::error expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, dd_edge &b)
{
  static const int nForests = 2;
  static forest* forests[nForests];
  forests[0] = a.getForest();
  forests[1] = b.getForest();
  op_info* opInfo = getOpInfo(op, forests, nForests);
  return opInfo == 0?
    compute_manager::UNKNOWN_OPERATION:
    apply(opInfo, a, b);
}


compute_manager::error expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, const dd_edge &b,
    dd_edge &c)
{
  static const int nForests = 3;
  static forest* forests[nForests];
  forests[0] = a.getForest();
  forests[1] = b.getForest();
  forests[2] = c.getForest();
  op_info* opInfo = getOpInfo(op, forests, nForests);
  return opInfo == 0?
    compute_manager::UNKNOWN_OPERATION:
    apply(opInfo, a, b, c);
}


compute_manager::error expert_compute_manager::setHashTablePolicy(
    bool chaining, unsigned size)
{
  if (size == 0) return compute_manager::TYPE_MISMATCH;
  if (cc->setPolicy(chaining, size))
    return compute_manager::SUCCESS;
  else
    return compute_manager::TYPE_MISMATCH;
}


void expert_compute_manager::showComputeTable(FILE* strm) const
{
  if (cc != 0) cc->show(strm);
}


int expert_compute_manager::getNumCacheEntries() const
{
  return (cc == 0)? 0: cc->getNumEntries();
}


void expert_compute_manager::removeStales()
{
  if (cc != 0) cc->removeStales();
}


void expert_compute_manager::clearComputeTable()
{
  if (cc != 0) cc->clear();
}


void expert_compute_manager::addBuiltinOp(const builtin_op_key& key,
  const operation* op, const forest* const* forests, int n)
{
  op_info entry(const_cast<operation*>(op), const_cast<forest**>(forests),
      n, cc);
  (*builtinOpEntries)[key] = entry;
#ifdef DEVELOPMENT_CODE
  std::map<builtin_op_key, op_info>::iterator itr = builtinOpEntries->find(key);
  if (itr == builtinOpEntries->end()) {
    // print the entries
    itr = builtinOpEntries->begin();
    while (itr != builtinOpEntries->end()) {
      printf("{");
      itr->first.print(stdout);
      itr->second.print(stdout);
      printf("}\n");
      ++itr;
    }
    assert(false);
  }
  assert(itr->second == entry);
#endif
}


op_info* expert_compute_manager::getOpInfo(compute_manager::op_code op,
    const forest* const* forests, int N)
{
  // search in built-in op entries
  builtin_op_key key(op, forests, N);
  std::map<builtin_op_key, op_info>::iterator itr = builtinOpEntries->find(key);
  if (itr != builtinOpEntries->end()) return &(itr->second);

  // add new built-in op entry
  if (N == 2) {
    // unary operations
    if (smart_cast<const expert_forest* const>(forests[0])->isMdd() &&
        smart_cast<const expert_forest* const>(forests[1])->isMtMdd() &&
        op == compute_manager::COPY) {
      // MDD to MTMDD
      // terminal true == 1, terminal false == 0
      addBuiltinOp(key, mdd_to_mtmdd::getInstance(), forests, N);
      return &(builtinOpEntries->find(key)->second);
    }
    if (smart_cast<const expert_forest* const>(forests[0])->isMtMdd() &&
        op == compute_manager::COPY) {
      if (smart_cast<const expert_forest* const>(forests[1])->isMdd()) {
        // MTMDD to MDD
        // terminal 0 == false, !0 == true
        addBuiltinOp(key, mtmdd_to_mdd::getInstance(), forests, N);
        return &(builtinOpEntries->find(key)->second);
      }
      if (forests[1]->getEdgeLabeling() == forest::EVPLUS ||
          forests[1]->getEdgeLabeling() == forest::EVTIMES) {
        // MTMDD to EVMDD (works for both EVPLUS and EVTIMES)
        // terminal 0 == false, !0 == true
        addBuiltinOp(key, mtmdd_to_evmdd::getInstance(), forests, N);
        return &(builtinOpEntries->find(key)->second);
      }
    }
    if (smart_cast<const expert_forest* const>(forests[0])->isMxd() &&
        smart_cast<const expert_forest* const>(forests[1])->isMtMxd() &&
        op == compute_manager::COPY) {
      // MXD to MTMXD
      // terminal true == 1, terminal false == 0
      addBuiltinOp(key, mxd_to_mtmxd::getInstance(), forests, N);
      return &(builtinOpEntries->find(key)->second);
    }
    if (smart_cast<const expert_forest* const>(forests[0])->isMtMxd() &&
        smart_cast<const expert_forest* const>(forests[1])->isMxd() &&
        op == compute_manager::COPY) {
      // MTMXD to MXD
      // terminal 0 == false, !0 == true
      addBuiltinOp(key, mtmxd_to_mxd::getInstance(), forests, N);
      return &(builtinOpEntries->find(key)->second);
    }
  }
  else if (N == 3) {

    // binary operations
    if (smart_cast<const expert_forest* const>(forests[0])->isMdd()) {

      if (smart_cast<const expert_forest* const>(forests[1])->isMdd() &&
          smart_cast<const expert_forest* const>(forests[2])->isMdd()) {
        // MDD binary operation
        switch (op) {
          case UNION:
            // Mdd union
            addBuiltinOp(key, mdd_union::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case INTERSECTION:
            // Mdd intersection
            addBuiltinOp(key, mdd_intersection::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case DIFFERENCE:
            // Mdd difference
            addBuiltinOp(key, mdd_difference::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          default:
            break;
        }
      }

      if (smart_cast<const expert_forest* const>(forests[1])->isMxd() &&
          smart_cast<const expert_forest* const>(forests[2])->isMdd()) {
        // MDD MXD image operations
        switch (op) {
          case POST_IMAGE:
            // Mdd Post-Image
            addBuiltinOp(key, mdd_post_image::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case PRE_IMAGE:
            // Mdd Pre-Image
            addBuiltinOp(key, mdd_pre_image::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case REACHABLE_STATES_DFS:
            // Mdd Reachable states using saturation-based algorithm
            assert(false);
          case REACHABLE_STATES_BFS:
            // Mdd Reachable states using traditional breadth-first algorithm
            addBuiltinOp(key, mdd_reachability_bfs::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          default:
            break;
        }
      }

    }
    else if (smart_cast<const expert_forest* const>(forests[0])->isMxd()) {

      if (smart_cast<const expert_forest* const>(forests[1])->isMxd() &&
          smart_cast<const expert_forest* const>(forests[2])->isMxd()) {
        // MXD binary operation
        switch (op) {
          case UNION:
            // Mxd union
            addBuiltinOp(key, mxd_union::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case INTERSECTION:
            // Mxd intersection
            addBuiltinOp(key, mxd_intersection::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case DIFFERENCE:
            // Mxd difference
            addBuiltinOp(key, mxd_difference::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          default:
            break;
        }
      }

    }
    else if (smart_cast<const expert_forest* const>(forests[0])->isMtMdd()) {

      if (smart_cast<const expert_forest* const>(forests[1])->isMtMdd()) {

        if (smart_cast<const expert_forest* const>(forests[2])->isMtMdd()) {
          // MTMDD binary operation
          switch (op) {
            case MIN:
              // MtMdd min
              addBuiltinOp(key, mtmdd_min::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case MAX:
              // MtMdd max
              addBuiltinOp(key, mtmdd_max::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case PLUS:
              // MtMdd plus
              addBuiltinOp(key, mtmdd_plus::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case MINUS:
              // MtMdd minus
              addBuiltinOp(key, mtmdd_minus::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case MULTIPLY:
              // MtMdd multiply
              addBuiltinOp(key, mtmdd_multiply::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case DIVIDE:
              // MtMdd divide
              addBuiltinOp(key, mtmdd_divide::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case EQUAL:
              // MtMdd ==
              addBuiltinOp(key, mtmdd_equal::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case NOT_EQUAL:
              // MtMdd !=
              addBuiltinOp(key, mtmdd_not_equal::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case LESS_THAN:
              // MtMdd <
              addBuiltinOp(key, mtmdd_less_than::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case LESS_THAN_EQUAL:
              // MtMdd <=
              addBuiltinOp(key, mtmdd_less_than_equal::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case GREATER_THAN:
              // MtMdd >
              addBuiltinOp(key, mtmdd_greater_than::getInstance(), forests, N);
              return &(builtinOpEntries->find(key)->second);
            case GREATER_THAN_EQUAL:
              // MtMdd >=
              addBuiltinOp(key, mtmdd_greater_than_equal::getInstance(),
                  forests, N);
              return &(builtinOpEntries->find(key)->second);
            default:
              break;
          }
        }
      } // MTMDD binary operations
      else if (smart_cast<const expert_forest* const>(forests[1])->isMtMxd()) {

        if (smart_cast<const expert_forest* const>(forests[2])->isMtMdd()) {
          // MTMDD-MTMXD image operation
          switch (op) {
          case POST_IMAGE:
            // MtMdd Post-Image
            addBuiltinOp(key, mtmdd_post_image::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case PRE_IMAGE:
            // MtMdd Pre-Image
            addBuiltinOp(key, mtmdd_pre_image::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case REACHABLE_STATES_DFS:
            // MtMdd Reachable states using saturation-based algorithm
            assert(false);
          case REACHABLE_STATES_BFS:
            // MtMdd Reachable states using traditional breadth-first algorithm
            assert(false);
          default:
            break;
          }
        } // MTMDD-MTMXD image operation

      }

    }
    else if (smart_cast<const expert_forest* const>(forests[0])->isMtMxd()) {

      if (smart_cast<const expert_forest* const>(forests[1])->isMtMxd() &&
          smart_cast<const expert_forest* const>(forests[2])->isMtMxd()) {
        // MTMXD binary operation
        switch (op) {
          case MIN:
            // MtMxd min
            addBuiltinOp(key, mtmxd_min::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case MAX:
            // MtMxd max
            addBuiltinOp(key, mtmxd_max::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case PLUS:
            // MtMxd plus
            addBuiltinOp(key, mtmxd_plus::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case MINUS:
            // MtMxd minus
            addBuiltinOp(key, mtmxd_minus::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case MULTIPLY:
            // MtMxd multiply
            addBuiltinOp(key, mtmxd_multiply::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case DIVIDE:
            // MtMxd divide
            addBuiltinOp(key, mtmxd_divide::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case EQUAL:
            // MtMxd ==
            addBuiltinOp(key, mtmxd_equal::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case NOT_EQUAL:
            // MtMxd !=
            addBuiltinOp(key, mtmxd_not_equal::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case LESS_THAN:
            // MtMxd <
            addBuiltinOp(key, mtmxd_less_than::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case LESS_THAN_EQUAL:
            // MtMxd <=
            addBuiltinOp(key, mtmxd_less_than_equal::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case GREATER_THAN:
            // MtMxd >
            addBuiltinOp(key, mtmxd_greater_than::getInstance(), forests, N);
            return &(builtinOpEntries->find(key)->second);
          case GREATER_THAN_EQUAL:
            // MtMxd >=
            addBuiltinOp(key, mtmxd_greater_than_equal::getInstance(),
                forests, N);
            return &(builtinOpEntries->find(key)->second);
          default:
            break;
        }
      }

    }

  }
  return 0;
}


op_info* expert_compute_manager::getOpInfo(const operation* op,
    const forest* const* forests, int N)
{
  // search in custom op entries
  custom_op_key key(op, forests, N);
  std::map<custom_op_key, op_info>::iterator itr = customOpEntries->find(key);
  if (itr == customOpEntries->end()) {
    // add new entry
    op_info entry(const_cast<operation*>(op), const_cast<forest**>(forests),
        N, cc);
    (*customOpEntries)[key] = entry;
    itr = customOpEntries->find(key);
#ifdef DEVELOPMENT_CODE
    assert (itr != customOpEntries->end());
    assert(itr->second == entry);
#endif
  }
  return &(itr->second);
}


