
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



#include "defines.h"
#include "operations/operation_ext.h"
#include "operations/reachability.h"
#include "operations/cardinality.h"
#include "operations/maxmin_range.h"
#include "operations/vect_matr.h"
#include "operations/cross.h"
#include "compute_cache.h"

#include "revision.h"

// ------------------ compute_manager --------------------------
MEDDLY::compute_manager::compute_manager() {}


MEDDLY::compute_manager::~compute_manager() {}


const char*
MEDDLY::compute_manager::getErrorCodeName(compute_manager::error e)
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


compute_manager* MEDDLY::getComputeManager()
{
  static expert_compute_manager* ecm = new expert_compute_manager();
  return ecm;
}


const char* MEDDLY::getLibraryInfo(int what)
{
  static char* title = 0;
  switch (what) {
    case 0:
      if (!title) {
        title = new char[80];
        if (REVISION_NUMBER) {
          snprintf(title, 80, 
            "%s version %s.%d (32-bit and 64-bit compatible)", 
            PACKAGE_NAME, VERSION, REVISION_NUMBER
          );
        } else {
          snprintf(title, 80, 
            "%s version %s (32-bit and 64-bit compatible)", 
            PACKAGE_NAME, VERSION
          );
        }
      }
      return title;

    case 1:
      return "Copyright (C) 2009, Iowa State University Research Foundation, Inc.";

    case 2:
      return "Released under the GNU Lesser General Public License, version 3";
 
    case 3:
      return "http://meddly.sourceforge.net/";

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

MEDDLY::expert_compute_manager::expert_compute_manager()
{
  // initialize compute cache
  cc = new compute_cache();
  // initialize builtinOpEntries and customOpEntries
  builtinOpEntries = new std::map<builtin_op_key, op_info>();
  customOpEntries = new std::map<custom_op_key, op_info>();
}


MEDDLY::expert_compute_manager::~expert_compute_manager()
{
  delete customOpEntries;
  delete builtinOpEntries;
  delete cc;
}


const char* MEDDLY::expert_compute_manager::getOperationName(
    compute_manager::op_code op) const
{
  switch(op) {
    case COPY:          return "Copy Edge";
    case CARDINALITY:   return "Cardinality";
    case UNION:         return "Union";
    case INTERSECTION:  return "Intersection";
    case DIFFERENCE:    return "Difference";
    case COMPLEMENT:    return "Complement";
    case MAX_RANGE:     return "Maximum range";
    case MIN_RANGE:     return "Minimum range";
    case PRE_IMAGE:     return "Pre-Image";
    case POST_IMAGE:    return "Post-Image";
    case REACHABLE_STATES_DFS:
                     return "Reachable States via Depth-First Search";
    case REACHABLE_STATES_BFS:
                     return "Reachable States via Breadth-First Search";
    case REVERSE_REACHABLE_DFS:
                     return "Reverse Reachable States via Depth-First Search";
    case REVERSE_REACHABLE_BFS:
                     return "Reverse Reachable States via Breadth-First Search";
    default: return "Unknown operation";
  }
}

template <class TYPE>
inline compute_manager::error
unary_apply(op_info* owner, const dd_edge &a, TYPE &b)
{
  // type check
  compute_manager::error err = owner->op->typeCheck(owner);
  return err != compute_manager::SUCCESS?
    err:
    owner->op->compute(owner, a, b);
}

compute_manager::error MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, dd_edge &b)
{
  return unary_apply(owner, a, b);
}

compute_manager::error MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, long &b)
{
  return unary_apply(owner, a, b);
}

compute_manager::error MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, double &b)
{
  return unary_apply(owner, a, b);
}

compute_manager::error MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, ct_object &b)
{
#ifdef HAVE_LIBGMP
  return unary_apply(owner, a, b);
#else
  return compute_manager::UNKNOWN_OPERATION;
#endif
}


compute_manager::error MEDDLY::expert_compute_manager::apply(op_info* owner,
    const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  // type check
  compute_manager::error err = owner->op->typeCheck(owner);
  return err != compute_manager::SUCCESS?
    err:
    owner->op->compute(owner, a, b, c);
}


template <class TYPE>
inline compute_manager::error 
unary_apply(expert_compute_manager* CM, compute_manager::op_code op, 
            const dd_edge &a, TYPE &b)
{
  static op_param plist[2];
  plist[0].set(a);
  plist[1].set(b);
  op_info* opInfo = CM->getOpInfo(op, plist, 2);
  return opInfo == 0?
    compute_manager::UNKNOWN_OPERATION:
    CM->apply(opInfo, a, b);
}

compute_manager::error MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, dd_edge &b)
{
  return unary_apply(this, op, a, b);
}

compute_manager::error MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, long &b)
{
  return unary_apply(this, op, a, b);
}

compute_manager::error MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, double &b)
{
  return unary_apply(this, op, a, b);
}

compute_manager::error MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, ct_object &b)
{
#ifdef HAVE_LIBGMP
  return unary_apply(this, op, a, b);
#else
  return compute_manager::UNKNOWN_OPERATION;
#endif
}


compute_manager::error MEDDLY::expert_compute_manager::apply(
    compute_manager::op_code op, const dd_edge &a, const dd_edge &b,
    dd_edge &c)
{
  static const int nForests = 3;
  static op_param plist[nForests];
  plist[0].set(a);
  plist[1].set(b);
  plist[2].set(c);
  op_info* opInfo = getOpInfo(op, plist, nForests);
  return opInfo == 0?
    compute_manager::UNKNOWN_OPERATION:
    apply(opInfo, a, b, c);
}


compute_manager::error MEDDLY::expert_compute_manager::setHashTablePolicy(
    bool chaining, unsigned size)
{
  if (size == 0) return compute_manager::TYPE_MISMATCH;
  if (cc->setPolicy(chaining, size))
    return compute_manager::SUCCESS;
  else
    return compute_manager::TYPE_MISMATCH;
}


void MEDDLY::expert_compute_manager::showComputeTable(FILE* strm) const
{
  if (cc != 0) cc->show(strm, true);
}


long MEDDLY::expert_compute_manager::getNumCacheEntries() const
{
  return (cc == 0)? 0: cc->getNumEntries();
}


void MEDDLY::expert_compute_manager::removeStales(op_info* op)
{
  if (cc != 0) cc->removeStales(op);
}


void MEDDLY::expert_compute_manager::clearComputeTable()
{
  if (cc != 0) cc->clear();
}


void MEDDLY::expert_compute_manager::addBuiltinOp(const builtin_op_key& key,
  const operation* op, const op_param* plist, int n)
{
  op_info entry(const_cast<operation*>(op), const_cast<op_param*>(plist),
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


op_info* MEDDLY::expert_compute_manager::getOpInfo(compute_manager::op_code op,
    const op_param* const plist, int N)
{
  // search in built-in op entries
  builtin_op_key key(op, plist, N);
  std::map<builtin_op_key, op_info>::iterator itr = builtinOpEntries->find(key);
  if (itr != builtinOpEntries->end()) return &(itr->second);

  const expert_forest* const f0 = plist[0].readForest();
  const expert_forest* const f1 = plist[1].readForest();
  const expert_forest* const f2 = plist[2].readForest();

  // add new built-in op entry
  if (N == 2) {
    // unary operations
    operation* opera = 0;
    switch (op) {
        case compute_manager::COMPLEMENT:
            if (f0->isMdd()) {
              if (f1->isMdd()) {
                opera = mdd_complement::getInstance();
              }
            } else if (f0->isMxd()) {
              if (f1->isMxd()) {
                opera = mxd_complement::getInstance();
              }
            }
            break;

        case compute_manager::CARDINALITY:
            opera = getCardinalityOperation(plist[0], plist[1]);
            break;

        case compute_manager::MAX_RANGE:
            opera = getMaxRangeOperation(plist[0], plist[1]);
            break;

        case compute_manager::MIN_RANGE:
            opera = getMinRangeOperation(plist[0], plist[1]);
            break;

        case compute_manager::CONVERT_TO_INDEX_SET:
            if (f0->isMdd() && f1->isEvplusMdd()) {
              opera = mdd_to_evplusmdd_index_set::getInstance();
            }
            break;

        case compute_manager::COPY:
            if (f0->isMdd()) {
              if (f1->isMtMdd()) {
                // MDD to MTMDD
                // terminal true == 1, terminal false == 0
                opera = mdd_to_mtmdd::getInstance();
              } // f1
            } else if (f0->isMtMdd()) {
              if (f1->isMdd()) {
                // MTMDD to MDD
                // terminal 0 == false, !0 == true
                opera = mtmdd_to_mdd::getInstance();
              } else if ( f1->getEdgeLabeling() == forest::EVPLUS ||
                          f1->getEdgeLabeling() == forest::EVTIMES) {
                // MTMDD to EVMDD (works for both EVPLUS and EVTIMES)
                // terminal 0 == false, !0 == true
                opera = mtmdd_to_evmdd::getInstance();
              } // f1
            } else if (f0->isMxd()) {
              if (f1->isMtMxd()) {
                // MXD to MTMXD
                // terminal true == 1, terminal false == 0
                opera = mxd_to_mtmxd::getInstance();
              } // f1
            } else if (f0->isMtMxd()) {
              if (f1->isMxd()) {
                // MTMXD to MXD
                // terminal 0 == false, !0 == true
                opera = mtmxd_to_mxd::getInstance();
              } // f1
            } // f0
            break;

        default:
            // not strictly necessary, but keeps compilers happy
            return 0;

    } // switch
    if (0==opera) return 0;
    addBuiltinOp(key, opera, plist, N);
    return &(builtinOpEntries->find(key)->second);

  } // End of Unary operations

  if (N == 3) {
    // binary operations
    operation* opera = 0;
    if (CROSS == op) {
      opera = getCrossOperation(plist[0]);
      if (0==opera) return 0;
      addBuiltinOp(key, opera, plist, N);
      return &(builtinOpEntries->find(key)->second);
    }

    if (f0->isMdd()) {
      if (f1->isMdd() &&
          f2->isMdd()) {
        // MDD binary operation
        switch (op) {
          case UNION:
            // Mdd union
            addBuiltinOp(key, mdd_union::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case INTERSECTION:
            // Mdd intersection
            addBuiltinOp(key, mdd_intersection::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case DIFFERENCE:
            // Mdd difference
            addBuiltinOp(key, mdd_difference::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          default:
            break;
        }
      }

      if (f1->isMxd() &&
          f2->isMdd()) {
        // MDD MXD image operations
        switch (op) {
          case POST_IMAGE:
            // Mdd Post-Image
            addBuiltinOp(key, mdd_post_image::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case PRE_IMAGE:
            // Mdd Pre-Image
            addBuiltinOp(key, mdd_pre_image::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case REACHABLE_STATES_DFS:
            // Mdd Reachable states using saturation-based algorithm
            addBuiltinOp(key, mdd_reachability_dfs::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case REACHABLE_STATES_BFS:
            // Mdd Reachable states using traditional breadth-first algorithm
            addBuiltinOp(key, mdd_reachability_bfs::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case REVERSE_REACHABLE_BFS:
            // Mdd Reverse reachable states using traditional
            // breadth-first algorithm
            addBuiltinOp(key, mdd_backward_reachability_bfs::getInstance(),
                plist, N);
            return &(builtinOpEntries->find(key)->second);
          case REVERSE_REACHABLE_DFS:
            // Mdd Reverse Reachable states using saturation-based algorithm
            addBuiltinOp(key, mdd_backward_reachability_dfs::getInstance(),
                plist, N);
            return &(builtinOpEntries->find(key)->second);
          default:
            break;
        }
      }
    }
    else if (f0->isMxd()) {

      if (f1->isMxd() &&
          f2->isMxd()) {
        // MXD binary operation
        switch (op) {
          case UNION:
            // Mxd union
            addBuiltinOp(key, mxd_union::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case INTERSECTION:
            // Mxd intersection
            addBuiltinOp(key, mxd_intersection::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case DIFFERENCE:
            // Mxd difference
            addBuiltinOp(key, mxd_difference::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          default:
            break;
        }
      }

    }
    else if (f0->isMtMdd()) {

      if (f1->isMtMdd()) {

        if (f2->isMtMdd()) {
          // MTMDD binary operation
          switch (op) {
            case MIN:
              // MtMdd min
              addBuiltinOp(key, mtmdd_min::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case MAX:
              // MtMdd max
              addBuiltinOp(key, mtmdd_max::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case PLUS:
              // MtMdd plus
              addBuiltinOp(key, mtmdd_plus::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case MINUS:
              // MtMdd minus
              addBuiltinOp(key, mtmdd_minus::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case MULTIPLY:
              // MtMdd multiply
              addBuiltinOp(key, mtmdd_multiply::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case DIVIDE:
              // MtMdd divide
              addBuiltinOp(key, mtmdd_divide::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case EQUAL:
              // MtMdd ==
              addBuiltinOp(key, mtmdd_equal::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case NOT_EQUAL:
              // MtMdd !=
              addBuiltinOp(key, mtmdd_not_equal::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case LESS_THAN:
              // MtMdd <
              addBuiltinOp(key, mtmdd_less_than::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case LESS_THAN_EQUAL:
              // MtMdd <=
              addBuiltinOp(key, mtmdd_less_than_equal::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case GREATER_THAN:
              // MtMdd >
              addBuiltinOp(key, mtmdd_greater_than::getInstance(), plist, N);
              return &(builtinOpEntries->find(key)->second);
            case GREATER_THAN_EQUAL:
              // MtMdd >=
              addBuiltinOp(key, mtmdd_greater_than_equal::getInstance(),
                  plist, N);
              return &(builtinOpEntries->find(key)->second);
            default:
              break;
          }
        }
      } // MTMDD binary operations
      else if (f1->isMtMxd()) {

        if (f2->isMtMdd()) {
          // MTMDD-MTMXD image operation
          switch (op) {
          case POST_IMAGE:
            // MtMdd Post-Image
            addBuiltinOp(key, mtmdd_post_image::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case PRE_IMAGE:
            // MtMdd Pre-Image
            addBuiltinOp(key, mtmdd_pre_image::getInstance(), plist, N);
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
    else if (f0->isMtMxd()) {

      if (f1->isMtMxd() &&
          f2->isMtMxd()) {
        // MTMXD binary operation
        switch (op) {
          case MIN:
            // MtMxd min
            addBuiltinOp(key, mtmxd_min::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MAX:
            // MtMxd max
            addBuiltinOp(key, mtmxd_max::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case PLUS:
            // MtMxd plus
            addBuiltinOp(key, mtmxd_plus::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MINUS:
            // MtMxd minus
            addBuiltinOp(key, mtmxd_minus::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MULTIPLY:
            // MtMxd multiply
            addBuiltinOp(key, mtmxd_multiply::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case DIVIDE:
            // MtMxd divide
            addBuiltinOp(key, mtmxd_divide::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case EQUAL:
            // MtMxd ==
            addBuiltinOp(key, mtmxd_equal::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case NOT_EQUAL:
            // MtMxd !=
            addBuiltinOp(key, mtmxd_not_equal::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case LESS_THAN:
            // MtMxd <
            addBuiltinOp(key, mtmxd_less_than::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case LESS_THAN_EQUAL:
            // MtMxd <=
            addBuiltinOp(key, mtmxd_less_than_equal::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case GREATER_THAN:
            // MtMxd >
            addBuiltinOp(key, mtmxd_greater_than::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case GREATER_THAN_EQUAL:
            // MtMxd >=
            addBuiltinOp(key, mtmxd_greater_than_equal::getInstance(),
                plist, N);
            return &(builtinOpEntries->find(key)->second);
          default:
            break;
        }
      }

    }
    else if (f0->isEvplusMdd()) {
      if (f1->isEvplusMdd() &&
          f2->isEvplusMdd()) {

        // EV+MDD binary operation
        switch (op) {
          case PLUS:
            // Ev+Mdd plus
            addBuiltinOp(key, evplusmdd_plus::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MINUS:
            // Ev+Mdd minus
            addBuiltinOp(key, evplusmdd_minus::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MULTIPLY:
            // Ev+Mdd multiply
            addBuiltinOp(key, evplusmdd_multiply::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
#if 0
          case DIVIDE:
            // Ev+Mdd divide
            addBuiltinOp(key, evplusmdd_divide::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MIN:
            // Ev+Mdd min
            addBuiltinOp(key, evplusmdd_min::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MAX:
            // Ev+Mdd max
            addBuiltinOp(key, evplusmdd_max::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case EQUAL:
            // Ev+Mdd ==
            addBuiltinOp(key, evplusmdd_equal::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case NOT_EQUAL:
            // Ev+Mdd !=
            addBuiltinOp(key, evplusmdd_not_equal::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case LESS_THAN:
            // Ev+Mdd <
            addBuiltinOp(key, evplusmdd_less_than::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case LESS_THAN_EQUAL:
            // Ev+Mdd <=
            addBuiltinOp(key, evplusmdd_less_than_equal::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case GREATER_THAN:
            // MtMdd >
            addBuiltinOp(key, evplusmdd_greater_than::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case GREATER_THAN_EQUAL:
            // MtMdd >=
            addBuiltinOp(key, evplusmdd_greater_than_equal::getInstance(),
                plist, N);
            return &(builtinOpEntries->find(key)->second);
#endif
          default:
            break;
        }
      }
    }
    else if (f0->isEvtimesMdd()) {
      if (f1->isEvtimesMdd() &&
          f2->isEvtimesMdd()) {

        // EV*MDD binary operation
        switch (op) {
          case PLUS:
            // Ev*Mdd plus
            addBuiltinOp(key, evtimesmdd_plus::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MINUS:
            // Ev*Mdd minus
            addBuiltinOp(key, evtimesmdd_minus::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MULTIPLY:
            // Ev*Mdd multiply
            addBuiltinOp(key, evtimesmdd_multiply::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
#if 0
          case DIVIDE:
            // Ev*Mdd divide
            addBuiltinOp(key, evtimesmdd_divide::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MIN:
            // Ev*Mdd min
            addBuiltinOp(key, evtimesmdd_min::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case MAX:
            // Ev*Mdd max
            addBuiltinOp(key, evtimesmdd_max::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
#endif
          case EQUAL:
            // Ev*Mdd ==
            addBuiltinOp(key, evtimesmdd_equal::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
#if 0
          case NOT_EQUAL:
            // Ev*Mdd !=
            addBuiltinOp(key, evtimesmdd_not_equal::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case LESS_THAN:
            // Ev*Mdd <
            addBuiltinOp(key, evtimesmdd_less_than::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case LESS_THAN_EQUAL:
            // Ev*Mdd <=
            addBuiltinOp(key, evtimesmdd_less_than_equal::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case GREATER_THAN:
            // Ev*Mdd >
            addBuiltinOp(key, evtimesmdd_greater_than::getInstance(), plist, N);
            return &(builtinOpEntries->find(key)->second);
          case GREATER_THAN_EQUAL:
            // Ev*Mdd >=
            addBuiltinOp(key, evtimesmdd_greater_than_equal::getInstance(),
                plist, N);
            return &(builtinOpEntries->find(key)->second);
#endif
          default:
            break;
        }
      }
    }
  }
  return 0;
}


op_info* MEDDLY::expert_compute_manager::getOpInfo(const operation* op,
    const op_param* const plist, int N)
{
  // search in custom op entries
  custom_op_key key(op, plist, N);
  std::map<custom_op_key, op_info>::iterator itr = customOpEntries->find(key);
  if (itr == customOpEntries->end()) {
    // add new entry
    op_info entry(const_cast<operation*>(op), 
        const_cast<op_param*>(plist), N, cc);
    (*customOpEntries)[key] = entry;
    itr = customOpEntries->find(key);
#ifdef DEVELOPMENT_CODE
    assert (itr != customOpEntries->end());
    assert(itr->second == entry);
#endif
  }
  return &(itr->second);
}

compute_manager::error 
MEDDLY::expert_compute_manager::vectorMatrixMultiply(double* y, const dd_edge &y_ind,
                      const double* x, const dd_edge &x_ind, const dd_edge &A)
{
  const expert_forest* const fy = (expert_forest*) y_ind.getForest();
  const expert_forest* const fx = (expert_forest*) x_ind.getForest();
  const expert_forest* const fA = (expert_forest*) A.getForest();

  if (
           (fy->getRangeType() != forest::INTEGER) 
        || (fy->isForRelations())
        || (fx->getRangeType() != forest::INTEGER)
        || (fx->isForRelations())
        || (fA->getRangeType() != forest::REAL)
        || (!fA->isForRelations())
      ) 
  {
    return TYPE_MISMATCH;
  }

  // A can't be fully reduced.
  if (forest::FULLY_REDUCED == fA->getReductionRule()) {
    return TYPE_MISMATCH;
  }

  // For now, fy and fx must be EV+MDDs.
  if (     (fy->getEdgeLabeling() != forest::EVPLUS) 
        || (fx->getEdgeLabeling() != forest::EVPLUS) )
  {
    return NOT_IMPLEMENTED;
  }

  //everyone must use the same domain
  if (      (fx->getDomain() != fy->getDomain()) 
        ||  (fx->getDomain() != fA->getDomain())  )
  {
    return TYPE_MISMATCH;
  }

  static op_param plist[5];
  plist[0].set(y);
  plist[1].set(y_ind);
  plist[2].set(x);
  plist[3].set(x_ind);
  plist[4].set(A);

  switch (fA->getEdgeLabeling()) {
    case forest::MULTI_TERMINAL:
      return vectorMatrixMult_evplus_mt(
        plist, fy->getDomain()->getNumVariables(), y, y_ind.getNode(),
        x, x_ind.getNode(), A.getNode()
      );

    case forest::EVTIMES:
      return vectorMatrixMult_evplus_evtimes(
        plist, fy->getDomain()->getNumVariables(), y, y_ind.getNode(),
        x, x_ind.getNode(), A.getNode()
      );

    default:
      return TYPE_MISMATCH;
  };


}


compute_manager::error 
MEDDLY::expert_compute_manager::matrixVectorMultiply(double* y, const dd_edge &y_ind,
                      const dd_edge &A, const double* x, const dd_edge &x_ind)
{
  const expert_forest* const fy = (expert_forest*) y_ind.getForest();
  const expert_forest* const fA = (expert_forest*) A.getForest();
  const expert_forest* const fx = (expert_forest*) x_ind.getForest();

  if (
           (fy->getRangeType() != forest::INTEGER) 
        || (fy->isForRelations())
        || (fx->getRangeType() != forest::INTEGER)
        || (fx->isForRelations())
        || (fA->getRangeType() != forest::REAL)
        || (!fA->isForRelations())
      ) 
  {
    return TYPE_MISMATCH;
  }

  // A can't be fully reduced.
  if (forest::FULLY_REDUCED == fA->getReductionRule()) {
    return TYPE_MISMATCH;
  }

  // For now, fy and fx must be EV+MDDs.
  if (     (fy->getEdgeLabeling() != forest::EVPLUS) 
        || (fx->getEdgeLabeling() != forest::EVPLUS) )
  {
    return NOT_IMPLEMENTED;
  }

  //everyone must use the same domain
  if (      (fx->getDomain() != fy->getDomain()) 
        ||  (fx->getDomain() != fA->getDomain())  )
  {
    return TYPE_MISMATCH;
  }

  static op_param plist[5];
  plist[0].set(y);
  plist[1].set(y_ind);
  plist[2].set(A);
  plist[3].set(x);
  plist[4].set(x_ind);

  switch (fA->getEdgeLabeling()) {
    case forest::MULTI_TERMINAL:
      return matrixVectorMult_evplus_mt(
        plist, fy->getDomain()->getNumVariables(), y, y_ind.getNode(),
        A.getNode(), x, x_ind.getNode()
      );

    case forest::EVTIMES:
      return matrixVectorMult_evplus_evtimes(
        plist, fy->getDomain()->getNumVariables(), y, y_ind.getNode(),
        A.getNode(), x, x_ind.getNode()
      );

    default:
      return TYPE_MISMATCH;
  };


}


