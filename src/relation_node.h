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

#ifndef MEDDLY_RELATION_NODE
#define MEDDLY_RELATION_NODE

namespace MEDDLY {
    class relation_node;
};

// ******************************************************************
// *                                                                *
// *                      relation_node  class                      *
// *                                                                *
// ******************************************************************

/** Pieces of an implicit relation.

 Each piece (this class) is a function of a single variable.
 The function specifies the "next state" for the given
 current state, but we are only allowed to depend on one
 state variable.

 This is an abstract base class.  The function is specified
 by deriving a class from this one, and specifying method
 nextOf().
 TBD - Need a bogus value for nextOf()...

 Additionally, you must specify method equals(), which is used
 to detect when two nodes actually represent the same function.

 */
class MEDDLY::relation_node {
public:
  /** Constructor.
   @param  signature   Hash for this node, such that
   two equal nodes must have the same
   signature.
   @param  level          Level affected.
   @param  down           Handle to a relation node below us.
   */
 relation_node(unsigned long signature, forest* f, int level, node_handle down, long e_val, long f_val, long inh);
 virtual ~relation_node();

  // the following should be inlined in meddly_expert.hh

   // Get the forest to which the node belongs
  expert_forest* getForest();

  /** A signature for this function.
   This helps class implicit_relation to detect duplicate
   functions (using a hash table where the signature
   is taken as the hash value).
   */
  unsigned long getSignature() const;

  /** The state variable affected by this part of the relation.
   */
  int getLevel() const;

  /** Pointer to the (ID of the) next piece of the relation.
   */
  rel_node_handle getDown() const;

  void setDown(rel_node_handle d);

  /** The unique ID for this piece.
   */
  rel_node_handle getID() const;

  /** Set the unique ID for this piece.
   */
  void setID(rel_node_handle ID);

  /** The token_update array for this piece.
   */
  long* getTokenUpdate() const;

  /** Set the token_update array for this piece.
   */
  void setTokenUpdate(long* token_update);

  /** Get the enable condition for this piece.
   */
  long getEnable() const;

  /** Set the enable condition for this piece.
   */
  void setEnable(long enable_val);

  /** Get the fire condition for this piece.
   */
  long getFire() const;

  /** Set the fire condition for this piece.
   */
  void setFire(long fire_val);

  /** Get the inhibit condition for this piece.
   */
  long getInhibit() const;

  /** Set the inhibit condition for this piece.
   */
  void setInhibit(long inh_val);

  /** The size of token_update array for this piece.
   */
  long getPieceSize() const;

  /** Set the size of token_update array for this piece.
   */
  void setPieceSize(long pS);

  /** Expand the tokenUpdate array as the variable increases
   */
  void expandTokenUpdate(long i);

  /** Set the tokenUpdate array at location i to val
   */
  void setTokenUpdateAtIndex(long i,long val);

  // the following must be provided in derived classes.

  /** If the variable at this level has value i,
   what should the new value be?
   */
  virtual long nextOf(long i);

  /** Determine if this node is equal to another one.
   */
  virtual bool equals(const relation_node* n) const;

private:
  unsigned long signature;
  expert_forest* f;
  int level;
  long enable;
  long fire;
  long inhibit;
  rel_node_handle down;
  rel_node_handle ID;
  long* token_update;
  long piece_size;

  // used by the hash table in implicit_relation
  relation_node* hash_chain;

  // friend class implicit_relation;
};  // class relation_node


// ******************************************************************
// *                                                                *
// *                  inlined relation_node methods                 *
// *                                                                *
// ******************************************************************


inline unsigned long
MEDDLY::relation_node::getSignature() const
{
  return signature;
}

inline MEDDLY::expert_forest*
MEDDLY::relation_node::getForest() {
  return f;
}

inline int
MEDDLY::relation_node::getLevel() const
{
  return level;
}

inline rel_node_handle
MEDDLY::relation_node::getDown() const
{
  return down;
}

inline void
MEDDLY::relation_node::setDown(rel_node_handle d)
{
  down = d;
}


inline rel_node_handle
MEDDLY::relation_node::getID() const
{
  return ID;
}

inline void
MEDDLY::relation_node::setID(rel_node_handle n_ID)
{
  ID=n_ID;
}

inline long
MEDDLY::relation_node::getFire() const
{
  return fire;
}

inline void
MEDDLY::relation_node::setFire(long fire_val)
{
  fire = fire_val;
}

inline long
MEDDLY::relation_node::getEnable() const
{
  return enable;
}

inline void
MEDDLY::relation_node::setEnable(long enable_val)
{
  enable = enable_val;
}

inline void
MEDDLY::relation_node::setInhibit(long inh)
{
  inhibit = inh;
}

inline long
MEDDLY::relation_node::getInhibit() const
{
  return inhibit;
}

inline long
MEDDLY::relation_node::getPieceSize() const
{
  return piece_size;
}

inline void
MEDDLY::relation_node::setPieceSize(long pS)
{
  piece_size=pS;
}

inline long*
MEDDLY::relation_node::getTokenUpdate() const
{
  return token_update;
}

inline
void
MEDDLY::relation_node::setTokenUpdate(long* n_token_update)
{
  token_update = n_token_update;
}



#endif
