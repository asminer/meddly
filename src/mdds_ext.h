
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



// TODO: Testing
// TODO: mdd_node_manager
// TODO: mxd_node_manager
// TODO: evmdd_node_manager
// TODO: mtmdd_node_manager
// TODO: mtmxd_node_manager (??)

/* 
  TODO: ensure this rule
  All extensions must over-ride either reduceNode() or normalizeAndReduceNode().
  Normalizing is only for edge-valued decision diagrams.
*/
  
#ifndef MDDS_EXTENSIONS_H
#define MDDS_EXTENSIONS_H


#include "mdds.h"

#define SORT_BUILD

class mdd_node_manager : public node_manager {
  public:

    mdd_node_manager(domain *d);
    ~mdd_node_manager();

    // Refer to meddly.h
    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    
    // Refer to meddly_expert.h
    virtual int reduceNode(int p);

  protected:

    // Create edge representing f(vlist[]) = term and store it in curr.
    void createEdge(const int* vlist, int term, dd_edge& e);

    // Create a node, at level k, whose ith index points to dptr.
    // If i is -1, all indices of the node will point to dptr.
    int createNode(int k, int i, int dptr);

#ifdef SORT_BUILD
    int sortBuild(int** list, int height, int begin, int end);
#endif

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);

    // Refer to meddly.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist, 
        const float* terms, int N, dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const;
};


class mtmdd_node_manager : public node_manager {
  public:

    mtmdd_node_manager(domain *d, forest::range_type t);
    ~mtmdd_node_manager();

    // Refer to meddly.h
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    
    // Refer to meddly_expert.h
    virtual int reduceNode(int p);

  protected:

    // Create edge representing f(vlist[]) = term and store it in curr.
    void createEdge(const int* vlist, int term, dd_edge& e);

    // Create a node, at level k, whose ith index points to dptr.
    // If i is -1, all indices of the node will point to dptr.
    int createNode(int k, int i, int dptr);

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);

    // Refer to meddly.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist, 
        const float* terms, int N, dd_edge &e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const;
};


class mxd_node_manager : public node_manager {
  // TODO: mxds can only be forest::IDENTITY_REDUCED
  public:

    mxd_node_manager(domain *d);
    ~mxd_node_manager();

    // Refer to meddly.h
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge& e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;

    // Refer to meddly_expert.h
    virtual int reduceNode(int p);

  protected:
  
    // Create edge representing vlist[] and store it in curr.
    void createEdge(const int* vlist, const int* vplist, dd_edge& e);

    // Create a node n, at level -k, whose jth index points to dptr.
    // Create a node m, at level +k, whose ith index points to n.
    // If i or j is -1, all indices of n will point to dptr and all of m will
    //    point to n.
    // If i or j is -2, simply returns dptr.
    int createNode(int k, int i, int j, int dptr);

    // Create a node, at level k, whose ith index points to dptr.
    // 0 <= i < level bound
    // Used by createNode(k, i, j, dptr)
    int createNode(int k, int i, int dptr);

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);

    // Refer to meddly.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist, 
        const float* terms, int N, dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const;
};


class mtmxd_node_manager : public node_manager {
  // TODO: mtmxds can only be forest::IDENTITY_REDUCED
  public:

    mtmxd_node_manager(domain *d, forest::range_type t);
    ~mtmxd_node_manager();

    // Refer to meddly.h
    error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge& e);
    error createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge& e);

    error createEdge(int val, dd_edge &e);
    error createEdge(float val, dd_edge &e);

    error evaluate(const dd_edge& f, const int* vlist, const int* vplist,
        int &term) const;
    error evaluate(const dd_edge& f, const int* vlist, const int* vplist,
        float &term) const;

    // Refer to meddly_expert.h
    virtual int reduceNode(int p);

  protected:

    // Creates an edge representing v[] vp[] = terminal node (not value),
    // and stores it in e.
    void createEdge(const int* v, const int* vp, int termNode, dd_edge& e);

    // Creates a top-level node representing {-1, -1, ..., -1} = terminal node
    // (not value), and returns it (returned node is already linked to.
    int createEdge(int termNode);

    // Create a node n, at level -k, whose jth index points to dptr.
    // Create a node m, at level +k, whose ith index points to n.
    // If i or j is -1, all indices of n will point to dptr and all of m will
    //    point to n.
    // If i or j is -2, simply returns dptr.
    int createNode(int k, int i, int j, int dptr);

    // Create a node, at level k, whose ith index points to dptr.
    // 0 <= i < level bound
    // Used by createNode(k, i, j, dptr)
    int createNode(int k, int i, int dptr);

#ifdef SORT_BUILD
    int sortBuild(int** list, int** plist, int* tlist,
        int height, int begin, int end);
#endif

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev);

    // Refer to meddly.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);

    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual error createEdge(bool val, dd_edge &e);

    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
};


class evmdd_node_manager : public node_manager {
  public:
    evmdd_node_manager(domain *d, forest::edge_labeling el);
    ~evmdd_node_manager();

    virtual int doOp(int a, int b) const = 0;
    virtual float doOp(float a, float b) const = 0;
    virtual int getDefaultEdgeValue() const = 0;
    virtual int getIdentityEdgeValue() const = 0;
    virtual void initEdgeValues(int p) = 0;

    // Refer to meddly.h
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e);
    virtual error createEdge(int val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const;
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const;

  protected:

    // Create edge representing vlist[], ending in a terminal value term
    // (or you could think of it as a vector with an edge-value equal to term),
    // and store it in curr.
    void createEdge(const int* vlist, int term, dd_edge &e);

    int createTempNode(int k, int sz, bool clear = true);

    // Create a node, at level k, whose ith index points to dptr with an
    // edge value ev. If i is -1, all indices of the node will point to dptr
    // with an edge value ev.
    //void createNode(int k, int i, int dptr, int ev, int& res, int& resEv);

    template <class T>
    void createNode(int k, int index, int dptr, T ev, int& res, T& resEv);

    // Create a sparse node, at level k, whose ith index points to dptr with
    // edge value ev.
    // 0 <= i < level bound.
    // Used by createNode(k, i, dptr, ev, res, resEv).
    //void createSparseNode(int k, int i, int dptr, int ev, int& res, int& resEv);

    template<class T>
    void createSparseNode(int k, int index, int dptr, T ev, int& res, T& resEv);

  public:

    // Refer to meddly_expert.h
    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual int reduceNode(int p);

    // Refer to meddly.h:
    // The following will either abort or return an error since they are not
    // applicable to this forest.

    virtual error createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual error createEdge(const int* const* vlist, const int* const* vplist, 
        const float* terms, int N, dd_edge &e);
    virtual error createEdge(bool val, dd_edge &e);
    virtual error createEdge(float val, dd_edge &e);
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const;
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const;
};


class evplusmdd_node_manager : public evmdd_node_manager {
  public:
    evplusmdd_node_manager(domain *d);
    ~evplusmdd_node_manager();

    virtual int doOp(int a, int b) const;
    virtual float doOp(float a, float b) const;
    virtual int getDefaultEdgeValue() const;
    virtual int getIdentityEdgeValue() const;
    virtual void initEdgeValues(int p);
    virtual void normalizeAndReduceNode(int& p, int& ev);
    virtual void normalizeAndReduceNode(int& p, float& ev) { }
};


class evtimesmdd_node_manager : public evmdd_node_manager {
  public:
    evtimesmdd_node_manager(domain *d);
    ~evtimesmdd_node_manager();
    
    virtual int doOp(int a, int b) const;
    virtual float doOp(float a, float b) const;
    virtual int getDefaultEdgeValue() const;
    virtual int getIdentityEdgeValue() const;
    virtual void initEdgeValues(int p);
    virtual void normalizeAndReduceNode(int& p, int& ev) { }
    virtual void normalizeAndReduceNode(int& p, float& ev);
};

























// ------------------------ Inline methods -----------------------------------

inline
int evplusmdd_node_manager::doOp(int a, int b) const
{
  return a + b;
}

inline
int evtimesmdd_node_manager::doOp(int a, int b) const
{
  return a * b;
}

inline
float evplusmdd_node_manager::doOp(float a, float b) const
{
  return a + b;
}

inline
float evtimesmdd_node_manager::doOp(float a, float b) const
{
  return a * b;
}

inline
int evplusmdd_node_manager::getDefaultEdgeValue() const
{
  return INF;
}

inline
int evtimesmdd_node_manager::getDefaultEdgeValue() const
{
  static int intNan = toInt(NAN);
  return intNan;
}

inline
int evtimesmdd_node_manager::getIdentityEdgeValue() const
{
  return toInt(1.0);
}

inline
int evplusmdd_node_manager::getIdentityEdgeValue() const
{
  return 0;
}


template<class T>
void evmdd_node_manager::createNode(int k, int index, int dptr, T ev,
    int& res, T& resEv)
{
  DCASSERT(dptr >= -1 && ev >= 0);
  DCASSERT(index >= -1);

  if (index == -1) {
    // if all edge values are the same for a evmdd node, it is
    // equivalent to say that they are all 0 (since we subtract the minimum
    // anyway).
    if (reductionRule == forest::FULLY_REDUCED) {
      res = sharedCopy(dptr);
      resEv = ev;
      return;
    }
    res = createTempNodeMaxSize(k, false);
    setAllDownPtrsWoUnlink(res, dptr);
    setAllEdgeValues(res, ev);
    normalizeAndReduceNode(res, resEv);
    return;
  }

  // a single downpointer points to dptr

#if 0

  res = createTempNode(k, index + 1);
  setDownPtrWoUnlink(res, index, dptr);
  setEdgeValue(res, index, ev);
  resEv = 0;
  normalizeAndReduceNode(res, resEv);

#else

  if (nodeStorage == FULL_STORAGE ||
      (nodeStorage == FULL_OR_SPARSE_STORAGE && index < 2)) {
    // Build a full node
    res = createTempNode(k, index + 1);
    setDownPtrWoUnlink(res, index, dptr);
    setEdgeValue(res, index, ev);
    normalizeAndReduceNode(res, resEv);
  }
  else {
    DCASSERT (nodeStorage == SPARSE_STORAGE ||
        (nodeStorage == FULL_OR_SPARSE_STORAGE && index >= 2));
    // Build a sparse node
    createSparseNode(k, index, dptr, ev, res, resEv);
    // res, resEv set by createSparseNode(..); do not call normalizeAndReduce
  }

#endif

}


template<class T>
void evmdd_node_manager::createSparseNode(int k, int index,
    int dptr, T ev, int& res, T& resEv)
{
  DCASSERT(k != 0);

  if (isTimeToGc()) { garbageCollect(); }

  DCASSERT(isValidLevel(k));
  CHECK_RANGE(0, index, getLevelSize(k));

  // get a location in address[] to store the node
  int p = getFreeNode(k);

#ifdef DEBUG_MDD_SET
  printf("%s: k: %d, index: %d, new p: %d\n", __func__, k, index, p);
  fflush(stdout);
#endif

  // fill in the location with p's address info
  address[p].level = k;
  address[p].offset = getHole(k, 7 /*4 + 3*/, true);
  address[p].cache_count = 0;

#ifdef DEBUG_MDD_SET
  printf("%s: offset: %d\n", __func__, address[p].offset);
  fflush(stdout);
#endif

  int* foo = level[mapLevel(k)].data + address[p].offset;
  foo[0] = 1;                     // #incoming
  foo[1] = getTempNodeId();
  foo[2] = -1;                    // size
  foo[3] = index;                 // index
  foo[4] = sharedCopy(dptr);      // downpointer
  foo[5] = getIdentityEdgeValue();// this is the only ev, set resEv = ev
  foo[6] = p;                     // pointer to this node in the address array

  resEv = ev;

#ifdef TRACK_DELETIONS
  cout << "Creating node " << p << "\n";
  cout.flush();
#endif

  incrNodesActivatedSinceGc();

  // search in unique table
  int q = find(p);
  if (getNull() == q) {
    // no duplicate found; insert into unique table
    insert(p);
    DCASSERT(getCacheCount(p) == 0);
    DCASSERT(find(p) == p);
    res = p;
  }
  else {
    // duplicate found; discard this node and return the duplicate
    unlinkNode(dptr);
    // code from deleteTempNode(p) adapted to work here
    {
      makeHole(k, getNodeOffset(p), 7);
      freeNode(p);
      if (level[mapLevel(k)].compactLevel) compactLevel(k);
    }
    res = sharedCopy(q);
  }
}



#endif
