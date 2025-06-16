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

#ifndef MEDDLY_RELFOREST
#define MEDDLY_RELFOREST

#include "defines.h"
#include "edge_value.h"
#include "domain.h"

namespace MEDDLY {
    class relforest;
};

// ******************************************************************
// *                                                                *
// *                        relforest  class                        *
// *                                                                *
// ******************************************************************

/**
    Interface for a generic relation, compatible with DD operations.
    Can be thought of as a forest of implicit matrix diagram nodes.

    Derived classes specify how the relation is stored internally,
    but the external interface requires:
        * each "node" has a unique identifer (ID).
        * zero ID always means the empty relation.
        * negative ID always means 'terminal node'

 */
class MEDDLY::relforest {
    public:
        /// An unpacked row or column of a relation node.
        struct slice {
            node_handle diag_down;
            edge_value  diag_ev;

            std::vector <unsigned>      index;
            std::vector <node_handle>   down;
            std::vector <edge_value>    values;
        };
    public:
        /// A generic relation node.
        class node {
            public:
                inline unsigned ID() const {
                    return my_ID;
                }
                inline unsigned level() const {
                    return my_level;
                }
                inline const relforest* getParent() const {
                    return parent;
                }

                /**
                    For the relation node, get the outgoing edges from i.
                    When viewed as a matrix, this obtains row i.
                    Default behavior is to throw a 'not implemented' error.
                        @param  i   (unprimed) variable value
                        @param  u   will be filled with outgoing edges

                        @return true, iff there is at least one outgoing edge
                */
                virtual bool outgoing(unsigned i, slice &u);


                /**
                    For the relation node, get the incoming edges to i.
                    When viewed as a matrix, this obtains column i.
                    Default behavior is to throw a 'not implemented' error.
                        @param  i   (primed) variable value
                        @param  u   will be filled with outgoing edges

                        @return true, iff there is at least one incoming edge
                */
                virtual bool incoming(unsigned i, slice &u);

            protected:
                node(const relforest* _parent);
                virtual ~node();

                inline void set_ID(unsigned id) {
                    my_ID = id;
                }
                inline void set_level(unsigned lvl) {
                    my_level = lvl;
                }

            private:
                const relforest* parent;
                unsigned my_ID;
                unsigned my_level;

        };
    public:
        relforest(domain *_D);
        virtual ~relforest();

        inline const domain* getDomain() const {
            return D;
        }
        inline domain* getDomain() {
            return D;
        }
        inline unsigned getNumVariables() const {
            return D->getNumVariables();
        }

        /**
            Build a relation node for the given handle.
                @param  ID      Node to return.
                                Undefined behavior if ID <= 0.
                @param  byrows  If true, we will never call 'incoming'
                                on the created node.
        */
        virtual node& getNode(node_handle ID, bool byrows = true) = 0;

        /**
            Recycle a relation node.
        */
        virtual void doneNode(node &) = 0;


        /// Unique relforest identifier; use 0 for 'no forest'
        inline unsigned FID() const { return fid; }

        /// Get the relforest with specified ID, or nullptr.
        static inline relforest* getForestWithID(unsigned id) {
            if (id >= all_forests.size()) return nullptr;
            return all_forests[id];
        }

    private:
        static inline void initStatics() {
            all_forests.clear();
            all_forests.push_back(nullptr);
        }
        static inline void freeStatics() {
            all_forests.clear();
        }
        static inline void registerForest(relforest* f) {
            ASSERT(__FILE__, __LINE__, f);
            f->fid = all_forests.size();
            all_forests.push_back(f);
            f->D->registerRelforest(f);
        }
        static inline void unregisterForest(relforest* f) {
            ASSERT(__FILE__, __LINE__, f);
            if (f->fid < all_forests.size()) {
                all_forests[f->fid] = nullptr;
            }
            f->D->unregisterRelforest(f);
        }

    private:
        unsigned fid;
        domain* D;

        static std::vector <relforest*> all_forests;
        friend class initializer_list;
};

#endif

