
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


/*! \file meddly.hh

    Implementation details for interface in meddly.h.
*/

#ifndef MEDDLY_HH
#define MEDDLY_HH

//---------------------- Inlines ---------------------------------------------

// MEDDLY::

inline MEDDLY::domain* MEDDLY::createDomain() { 
  return createDomain((variable**) 0, 0);
}

#ifdef __GMP_H__
inline void MEDDLY::apply(const unary_opname* op, const dd_edge &a, mpz_t &c) {
  ct_object& x = get_mpz_wrapper();
  apply(op, a, HUGEINT, x);
  unwrap(x, c);
}
#endif

// error::

inline MEDDLY::error::code MEDDLY::error::getCode() const { return errcode; }
inline const char* MEDDLY::error::getFile() const { return fname; }
inline int MEDDLY::error::getLine() const { return lineno; }

// ******************************************************************
// *                                                                *
// *                         memstats class                         *
// *                                                                *
// ******************************************************************

inline void MEDDLY::memstats::incMemUsed(size_t b)
{
// #ifdef TRACK_LOCAL_MEMORY
  memory_used += b;
  if (memory_used > peak_memory_used) {
    peak_memory_used = memory_used;
  }
// #endif
// #ifdef TRACK_GLOBAL_MEMORY
  global_memory_used += b;
  if (global_memory_used > global_peak_used) {
    global_peak_used = global_memory_used;
  }
// #endif
}

inline void MEDDLY::memstats::decMemUsed(size_t b)
{
// #ifdef TRACK_LOCAL_MEMORY
  MEDDLY_DCASSERT(memory_used >= b);
  memory_used -= b;
// #endif
// #ifdef TRACK_GLOBAL_MEMORY
  MEDDLY_DCASSERT(global_memory_used >= b);
  global_memory_used -= b;
// #endif
}

inline void MEDDLY::memstats::incMemAlloc(size_t b)
{
// #ifdef TRACK_LOCAL_MEMORY
  memory_alloc += b;
  if (memory_alloc > peak_memory_alloc) {
    peak_memory_alloc = memory_alloc;
  }
// #endif
// #ifdef TRACK_GLOBAL_MEMORY
  global_memory_alloc += b;
  if (global_memory_alloc > global_peak_alloc) {
    global_peak_alloc = global_memory_alloc;
  }
// #endif
}

inline void MEDDLY::memstats::decMemAlloc(size_t b)
{
// #ifdef TRACK_LOCAL_MEMORY
  MEDDLY_DCASSERT(memory_alloc >= b);
  memory_alloc -= b;
// #endif
// #ifdef TRACK_GLOBAL_MEMORY
  MEDDLY_DCASSERT(global_memory_alloc >= b);
  global_memory_alloc -= b;
// #endif
}

inline void MEDDLY::memstats::zeroMemUsed()
{
  MEDDLY_DCASSERT(global_memory_used >= memory_used);
  global_memory_used -= memory_used;
  memory_used = 0;
}

inline void MEDDLY::memstats::zeroMemAlloc()
{
  MEDDLY_DCASSERT(global_memory_alloc >= memory_alloc);
  global_memory_alloc -= memory_alloc;
  memory_alloc = 0;
}

inline size_t MEDDLY::memstats::getMemUsed() const
{
  return memory_used;
}

inline size_t MEDDLY::memstats::getMemAlloc() const
{
  return memory_alloc;
}

inline size_t MEDDLY::memstats::getPeakMemUsed() const
{
  return peak_memory_used;
}

inline size_t MEDDLY::memstats::getPeakMemAlloc() const
{
  return peak_memory_alloc;
}


inline size_t MEDDLY::memstats::getGlobalMemUsed()
{
  return global_memory_used;
}

inline size_t MEDDLY::memstats::getGlobalMemAlloc()
{
  return global_memory_alloc;
}

inline size_t MEDDLY::memstats::getGlobalPeakMemUsed()
{
  return global_peak_used;
}

inline size_t MEDDLY::memstats::getGlobalPeakMemAlloc()
{
  return global_peak_alloc;
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                          forest class                          *
// *                                                                *
// *                                                                *
// ******************************************************************


inline void MEDDLY::forest::policies::setFullStorage() { 
  storage_flags = ALLOW_FULL_STORAGE; 
}

inline void MEDDLY::forest::policies::setSparseStorage() { 
  storage_flags = ALLOW_SPARSE_STORAGE;
}

inline void MEDDLY::forest::policies::setFullOrSparse() { 
  storage_flags = ALLOW_FULL_STORAGE | ALLOW_SPARSE_STORAGE;
}

inline void MEDDLY::forest::policies::setFullyReduced() {
  reduction = FULLY_REDUCED;
}

inline void MEDDLY::forest::policies::setQuasiReduced() {
  reduction = QUASI_REDUCED;
}

inline void MEDDLY::forest::policies::setIdentityReduced() {
  reduction = IDENTITY_REDUCED;
}

inline void MEDDLY::forest::policies::setZeroSuppressionReduced() {
  reduction = ZERO_SUPPRESSION_REDUCED;
}

inline void MEDDLY::forest::policies::setNeverDelete() {
  deletion = NEVER_DELETE;
}

inline void MEDDLY::forest::policies::setOptimistic() {
  deletion = OPTIMISTIC_DELETION;
}

inline void MEDDLY::forest::policies::setPessimistic() {
  deletion = PESSIMISTIC_DELETION;
}

inline void MEDDLY::forest::policies::setLowestInversion() {
  reorder = reordering_type::LOWEST_INVERSION;
}

inline void MEDDLY::forest::policies::setHighestInversion() {
  reorder = reordering_type::HIGHEST_INVERSION;
}

inline void MEDDLY::forest::policies::setSinkDown() {
  reorder = reordering_type::SINK_DOWN;
}

inline void MEDDLY::forest::policies::setBringUp() {
  reorder = reordering_type::BRING_UP;
}

inline void MEDDLY::forest::policies::setLowestCost() {
  reorder = reordering_type::LOWEST_COST;
}

inline void MEDDLY::forest::policies::setLowestMemory() {
  reorder = reordering_type::LOWEST_MEMORY;
}

inline void MEDDLY::forest::policies::setRandom() {
  reorder = reordering_type::RANDOM;
}

inline void MEDDLY::forest::policies::setLARC() {
  reorder = reordering_type::LARC;
}

inline void MEDDLY::forest::policies::setVarSwap() {
  swap = variable_swap_type::VAR;
}

inline void MEDDLY::forest::policies::setLevelSwap() {
  swap = variable_swap_type::LEVEL;
}

// end of struct policies

// forest::statset::
inline void MEDDLY::forest::statset::incActive(long b) {
  active_nodes += b;
  if (active_nodes > peak_active) 
    peak_active = active_nodes;
  MEDDLY_DCASSERT(active_nodes >= 0);
}
inline void MEDDLY::forest::statset::decActive(long b) {
  active_nodes -= b;
  MEDDLY_DCASSERT(active_nodes >= 0);
}
/*
inline void MEDDLY::forest::statset::incMemUsed(long b) {
  memory_used += b;
  if (memory_used > peak_memory_used) 
    peak_memory_used = memory_used;
  MEDDLY_DCASSERT(memory_used >= 0);
}
inline void MEDDLY::forest::statset::decMemUsed(long b) {
  memory_used -= b;
  MEDDLY_DCASSERT(memory_used >= 0);
}
inline void MEDDLY::forest::statset::incMemAlloc(long b) {
  memory_alloc += b;
  if (memory_alloc > peak_memory_alloc) 
    peak_memory_alloc = memory_alloc;
  MEDDLY_DCASSERT(memory_alloc >= 0);
}
inline void MEDDLY::forest::statset::decMemAlloc(long b) {
  memory_alloc -= b;
  MEDDLY_DCASSERT(memory_alloc >= 0);
}
inline void MEDDLY::forest::statset::sawUTchain(int c) {
  if (c > max_UT_chain)
    max_UT_chain = c;
}
*/
// end of forest::statset

// forest::logger::
inline bool MEDDLY::forest::logger::recordingNodeCounts() const {
  return node_counts;
}

inline void MEDDLY::forest::logger::recordNodeCounts() {
  if (nfix) node_counts = true;
}

inline void MEDDLY::forest::logger::ignoreNodeCounts() {
  if (nfix) node_counts = false;
}

inline bool MEDDLY::forest::logger::recordingTimeStamps() const {
  return time_stamps;
}

inline void MEDDLY::forest::logger::recordTimeStamps() {
  if (nfix) time_stamps = true;
}

inline void MEDDLY::forest::logger::ignoreTimeStamps() {
  if (nfix) time_stamps = false;
}

inline void MEDDLY::forest::logger::fixLogger() {
  nfix = false;
}

// end of forest::logger

// forest::

inline int MEDDLY::forest::downLevel(int k) {
  return (k>0) ? (-k) : (-k-1);
}
inline int MEDDLY::forest::upLevel(int k) {
  return (k<0) ? (-k) : (-k-1);
}

inline const MEDDLY::forest::policies& MEDDLY::forest::getDefaultPoliciesMDDs()
{
  return mddDefaults;
}

inline const MEDDLY::forest::policies& MEDDLY::forest::getDefaultPoliciesMXDs()
{
  return mxdDefaults;
}

inline void MEDDLY::forest::setDefaultPoliciesMDDs(const policies& p)
{
  mddDefaults = p;
}

inline void MEDDLY::forest::setDefaultPoliciesMXDs(const policies& p)
{
  mxdDefaults = p;
}

inline unsigned MEDDLY::forest::FID() const {
  return fid;
}

inline const MEDDLY::domain* MEDDLY::forest::getDomain() const {
  return d;
}

inline MEDDLY::domain* MEDDLY::forest::useDomain() {
  return d;
}

inline bool MEDDLY::forest::isForRelations() const {
  return isRelation;
}

inline MEDDLY::forest::range_type MEDDLY::forest::getRangeType() const {
  return rangeType;
}

inline bool MEDDLY::forest::isRangeType(MEDDLY::forest::range_type r) const {
  return r == rangeType;
}

inline MEDDLY::forest::edge_labeling MEDDLY::forest::getEdgeLabeling() const {
  return edgeLabel;
}

inline bool MEDDLY::forest::isMultiTerminal() const {
  return MULTI_TERMINAL == edgeLabel;
}

inline bool MEDDLY::forest::isEVPlus() const {
  return EVPLUS == edgeLabel;
}

inline bool MEDDLY::forest::isIndexSet() const {
  return INDEX_SET == edgeLabel;
}

inline bool MEDDLY::forest::isEVTimes() const {
  return EVTIMES == edgeLabel;
}

inline bool MEDDLY::forest::matches(bool isR, MEDDLY::forest::range_type rt,
  MEDDLY::forest::edge_labeling el) const {
  return (isRelation == isR) && (rangeType == rt) && (edgeLabel == el);
}

inline const MEDDLY::forest::policies& MEDDLY::forest::getPolicies() const {
  return deflt;
}

inline MEDDLY::forest::policies& MEDDLY::forest::getPolicies() {
  return deflt;
}

inline MEDDLY::forest::policies::reduction_rule MEDDLY::forest::getReductionRule() const {
  return deflt.reduction;
}


inline bool MEDDLY::forest::isFullyReduced() const {
    return MEDDLY::forest::policies::FULLY_REDUCED == deflt.reduction;
}

inline bool MEDDLY::forest::isQuasiReduced() const {
    return MEDDLY::forest::policies::QUASI_REDUCED == deflt.reduction;
}

inline bool MEDDLY::forest::isIdentityReduced() const {
    return MEDDLY::forest::policies::IDENTITY_REDUCED == deflt.reduction;
}

inline bool MEDDLY::forest::isZeroSuppressionReduced() const {
    return MEDDLY::forest::policies::ZERO_SUPPRESSION_REDUCED == deflt.reduction;
}

inline bool MEDDLY::forest::isUserDefinedReduced() const {
    return MEDDLY::forest::policies::USER_DEFINED == deflt.reduction;
}

inline bool MEDDLY::forest::isEdgeSpecificReduced() const {
    return MEDDLY::forest::policies::EDGE_SPECIFIC == deflt.reduction;
}

inline int* MEDDLY::forest::getLevelReductionRule() const{
    return level_reduction_rule;
    
}

inline bool MEDDLY::forest::isFullyReduced(int k) const {
    if(k<0)
        return getLevelReductionRule()[2*(-k)]==-1;
    else
        return getLevelReductionRule()[2*(k) - 1]==-1;
    
}

inline bool MEDDLY::forest::isQuasiReduced(int k) const {
    if(k<0)
        return getLevelReductionRule()[2*(-k)]==-2;
    else
        return getLevelReductionRule()[2*(k) - 1]==-2;
}

inline bool MEDDLY::forest::isIdentityReduced(int k) const {
    if(k<0)
        return getLevelReductionRule()[2*(-k)]==-3;
    else
        return false;
}



inline MEDDLY::node_storage_flags MEDDLY::forest::getNodeStorage() const {
  return deflt.storage_flags;
}

inline MEDDLY::forest::policies::node_deletion MEDDLY::forest::getNodeDeletion() const {
  return deflt.deletion;
}

inline bool MEDDLY::forest::isPessimistic() const {
  return MEDDLY::forest::policies::PESSIMISTIC_DELETION == deflt.deletion;
}

inline bool MEDDLY::forest::areSparseNodesEnabled() const {
  return MEDDLY::forest::policies::ALLOW_SPARSE_STORAGE & deflt.storage_flags;
}

inline bool MEDDLY::forest::areFullNodesEnabled() const {
  return MEDDLY::forest::policies::ALLOW_FULL_STORAGE & deflt.storage_flags;
}

inline const MEDDLY::forest::statset& MEDDLY::forest::getStats() const { 
  return stats; 
}

inline const MEDDLY::memstats& MEDDLY::forest::getMemoryStats() const {
  return mstats;
}

inline long MEDDLY::forest::getCurrentNumNodes() const {
  return stats.active_nodes;
}

inline size_t MEDDLY::forest::getCurrentMemoryUsed() const {
  return mstats.getMemUsed();
}

inline size_t MEDDLY::forest::getCurrentMemoryAllocated() const {
  return mstats.getMemAlloc();
}

inline long MEDDLY::forest::getPeakNumNodes() const {
  return stats.peak_active;
}

inline size_t MEDDLY::forest::getPeakMemoryUsed() const {
  return mstats.getPeakMemUsed();
}

inline size_t MEDDLY::forest::getPeakMemoryAllocated() const {
  return mstats.getPeakMemAlloc();
}

inline bool MEDDLY::forest::isMarkedForDeletion() const {
  return is_marked_for_deletion;
}

inline void MEDDLY::forest::createEdgeForVar(int vh, bool pr, dd_edge& a) {
  switch (rangeType) {
    case BOOLEAN:   createEdgeForVar(vh, pr, (bool*)  0, a);  break;
    case INTEGER:   createEdgeForVar(vh, pr, (long*)  0, a);  break;
    case REAL:      createEdgeForVar(vh, pr, (float*) 0, a);  break;
    default:        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
  }
};

inline void MEDDLY::forest::setLogger(logger* L, const char* name) {
  theLogger = L;
  if (theLogger) theLogger->logForestInfo(this, name);
}

// forest::edge_visitor::
inline void MEDDLY::forest::visitRegisteredEdges(edge_visitor &ev) {
  for (unsigned i = 0; i < firstFree; ++i) {
    if (edge[i].edge) ev.visit(*(edge[i].edge));
  }
}

// end of class forest


// ******************************************************************
// *                                                                *
// *                                                                *
// *                         variable class                         *
// *                                                                *
// *                                                                *
// ******************************************************************

inline int MEDDLY::variable::getBound(bool primed) const { 
  return primed ? pr_bound : un_bound; 
}
inline const char* MEDDLY::variable::getName() const { return name; }
inline bool MEDDLY::variable::isExtensible() const { 
  return is_extensible;
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                         variable order class                   *
// *                                                                *
// *                                                                *
// ******************************************************************

inline int MEDDLY::variable_order::getVarByLevel(int level) const 
{
  MEDDLY_DCASSERT(level>=0);
  return level2var[size_t(level)];
}

inline int MEDDLY::variable_order::getLevelByVar(int var) const 
{
  MEDDLY_DCASSERT(var>=0);
  return var2level[size_t(var)];
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                          domain class                          *
// *                                                                *
// *                                                                *
// ******************************************************************

inline int MEDDLY::domain::getNumVariables() const { return nVars; }

inline int
MEDDLY::domain::getVariableBound(int lev, bool prime) const {
  return vars[lev]->getBound(prime);
}

inline const MEDDLY::variable* MEDDLY::domain::getVar(int var) const {
  return vars[var];
}

inline MEDDLY::variable* MEDDLY::domain::useVar(int lev) { return vars[lev]; }

inline bool MEDDLY::domain::hasForests() const { return forests; }

inline bool MEDDLY::domain::isMarkedForDeletion() const {
  return is_marked_for_deletion;
}

inline std::shared_ptr<const MEDDLY::variable_order> MEDDLY::domain::makeDefaultVariableOrder() {
  return default_var_order;
}

inline int MEDDLY::domain::ID() const { return my_index; }


// ******************************************************************
// *                                                                *
// *                                                                *
// *                         dd_edge  class                         *
// *                                                                *
// *                                                                *
// ******************************************************************


inline void MEDDLY::dd_edge::clear() {
  assert(index != -1);
  set(0);
  raw_value = 0;
}

inline MEDDLY::forest* MEDDLY::dd_edge::getForest() const { return parent; }

inline MEDDLY::node_handle MEDDLY::dd_edge::getNode() const { return node; }

inline double MEDDLY::dd_edge::getCardinality() const {
  double c;
  apply(CARDINALITY, *this, c);
  return c;
}

inline bool MEDDLY::dd_edge::operator==(const MEDDLY::dd_edge& e) const {
  if (parent != e.parent) return false;
  return (node == e.node) && (raw_value == e.raw_value);
}

inline bool MEDDLY::dd_edge::operator!=(const MEDDLY::dd_edge& e) const {
  return !(*this == e);
}

inline const MEDDLY::dd_edge
MEDDLY::dd_edge::operator+(const MEDDLY::dd_edge& e) const {
  return dd_edge(*this) += e;
}

inline const MEDDLY::dd_edge MEDDLY::dd_edge::operator*(const MEDDLY::dd_edge& e) const {
  return dd_edge(*this) *= e;
}

inline const MEDDLY::dd_edge MEDDLY::dd_edge::operator-(const MEDDLY::dd_edge& e) const {
  return dd_edge(*this) -= e;
}

inline void MEDDLY::dd_edge::setIndex(int ind) { index = ind; }
inline int MEDDLY::dd_edge::getIndex() const { return index; }

inline void MEDDLY::dd_edge::orphan() {
  parent = 0;
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                        enumerator class                        *
// *                                                                *
// *                                                                *
// ******************************************************************


// enumerator::iterator::
inline bool MEDDLY::enumerator::iterator::build_error() const {
  return 0==F;
}

inline int MEDDLY::enumerator::iterator::levelChanged() const {
  return level_change;
}

inline const int* MEDDLY::enumerator::iterator::getAssignments() const {
  return index;
}

// enumerator::
inline MEDDLY::enumerator::operator bool() const {
  return is_valid;
}

inline void MEDDLY::enumerator::operator++() {
  is_valid &= I->next();
}

inline const int* MEDDLY::enumerator::getAssignments() const {
  if (I && is_valid) return I->getAssignments(); else return 0;
}

inline const int* MEDDLY::enumerator::getPrimedAssignments() const {
  if (I && is_valid) return I->getPrimedAssignments(); else return 0;
}

inline void MEDDLY::enumerator::getValue(int &v) const {
  if (I && is_valid) I->getValue(v);
}

inline void MEDDLY::enumerator::getValue(long &v) const {
  if (I && is_valid) I->getValue(v);
}

inline void MEDDLY::enumerator::getValue(float &v) const {
  if (I && is_valid) I->getValue(v);
}

inline int MEDDLY::enumerator::levelChanged() const {
  if (I) return I->levelChanged();
  return 0;
}

inline MEDDLY::enumerator::type MEDDLY::enumerator::getType() const {
  return T;
}

#endif
