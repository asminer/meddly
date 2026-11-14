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

#include <iostream>
#include <cassert>

#include "meddly/meddly.h"
#include "utils/rsgen_model.h"

// #define CHECK_SWAP_EVENT

// ******************************************************************

class swapmodel : public rsgen_model {
        template <bool FESWN>   // for each slot which number?
        class swap_adjacent : public rsgen_model::event {
                int num_slots;  // total number of slots
                int top_slot;   // larger slot number
            public:
                /**
                    An event to exchange numbers in adjacent slots.
                        @param N    Number of slots
                        @param k    Larger slot.
                */
                swap_adjacent(int N, int k);

                virtual bool fire(const MEDDLY::minterm &curr,
                        MEDDLY::minterm &);

                virtual bool build_nsf(MEDDLY::dd_edge &R);
        };
    public:
        swapmodel();

        virtual bool process_arg(unsigned &i, const char** argv);
        virtual void show_usage(std::ostream &out, const char* arg0);

        virtual bool build_model(MEDDLY::output &out);
        virtual void build_initial(MEDDLY::minterm &init);

        inline bool has_dimensions() const {
            return N>0;
        }

    protected:
        virtual bool known_cardinality(double &card) const;
#ifdef HAVE_LIBGMP
        virtual bool known_cardinality(mpz_t card) const;
#endif

    private:
        int N;
        bool for_each_slot_which_number;
};

// ******************************************************************
// *                swapmodel::swap_adjacent methods                *
// ******************************************************************

template <bool FESWN>
swapmodel::swap_adjacent<FESWN>::swap_adjacent(int N, int k)
{
    assert(N>0);
    assert(k >= 1);
    assert(k < N);

    num_slots = N;
    if (FESWN) {
        top_slot = k+1;
    } else {
        top_slot = k;
    }
}

// ******************************************************************

template <bool FESWN>
bool swapmodel::swap_adjacent<FESWN>::fire(const MEDDLY::minterm &curr,
        MEDDLY::minterm &next)
{
    if (FESWN) {
        /*
            For each slot which number.
            Exchange top_slot with top_slot-1
        */
        for (unsigned i=next.getNumVars(); i; i--) {
            next.setVar(i, curr.getVar(i));
        }
        next.setVar(top_slot, curr.getVar(top_slot-1));
        next.setVar(top_slot-1, curr.getVar(top_slot));
    } else {
        /*
            For each number, which slot is it in.
            Check all numbers. If it is in slot top_slot, decrement.
            If it is in slot top_slot-1, increment. Otherwise no change.
        */
        for (unsigned i=next.getNumVars(); i; i--) {
            if (curr.getVar(i) == top_slot) {
                next.setVar(i, top_slot-1);
                continue;
            }
            if (curr.getVar(i) == top_slot-1) {
                next.setVar(i, top_slot);
                continue;
            }
            next.setVar(i, curr.getVar(i));
        }
    }
    return true;
}

// *****************************************************************

template <bool FESWN>
bool swapmodel::swap_adjacent<FESWN>::build_nsf(MEDDLY::dd_edge &R)
{
    using namespace MEDDLY;

    forest* F = R.getForest();

    edge_value v;
    node_handle p;

    if (FESWN) {
        /*
            All levels are identity, except top_slot and top_slot-1,
            and their primed counterparts.
            There we need to exchange values.

            We can build this bottom-up in one go:
                for i
                    (x_top_slot == i)
                    for j
                        x'_top_slot = j
                        x_top_slot-1 = j
                        x'_top_slot-1 = i
        */
        node_handle one = F->makeIdentitiesTo(
            F->handleForValue(1), 0, top_slot-2, ~0
        );

        //
        // node at level x_t
        //
        unpacked_node* Utu = unpacked_node::newWritable(F,
                top_slot, num_slots, FULL_ONLY);
        for (int i=0; i<num_slots; i++) {
            //
            // node at level x'_t
            //
            unpacked_node* Utp = unpacked_node::newWritable(F,
                    -top_slot, num_slots, FULL_ONLY);
            for (int j=0; j<num_slots; j++) {
                if (j==i) {
                    Utp->setFull(j, 0);
                    continue;
                }
                //
                // singleton node at level x'_s  (s = t-1)
                //
                unpacked_node* Usp = unpacked_node::newWritable(F,
                        -(top_slot-1), 1, SPARSE_ONLY);
                Usp->setSparse(0, i, F->linkNode(one));
                F->createReducedNode(Usp, v, p);
                assert(v.isVoid());

                //
                // singleton node at level x_s
                //
                unpacked_node* Usu = unpacked_node::newWritable(F,
                        top_slot-1, 1, SPARSE_ONLY);
                Usu->setSparse(0, j, p);
                F->createReducedNode(Usu, v, p);
                assert(v.isVoid());

                // attach to x'_t node
                Utp->setFull(j, p);
            }
            F->createReducedNode(Utp, v, p);
            assert(v.isVoid());
            Utu->setFull(i, p);
        } // for i

        F->createReducedNode(Utu, v, p);
        assert(v.isVoid());
        p = F->makeIdentitiesTo(p, top_slot, num_slots, ~0);

        R.set(p);
        F->unlinkNode(one);
        return true;
    } else {
        /*
            All levels are the same:
                value top_slot becomes top_slot-1
                value top_slot-1 becomes top_slot
                all other values remain unchanged
        */
        node_handle down = F->handleForValue(1);
        F->linkNode(down);

        for (int k=1; k<=num_slots; k++) {
            unpacked_node* nk = unpacked_node::newWritable(F,
                    k, num_slots, FULL_ONLY);

            for (int i=0; i<num_slots; i++) {
                if (i == top_slot) {
                    unpacked_node* nkp = unpacked_node::newWritable(F,
                            -k, 1, SPARSE_ONLY);
                    nkp->setSparse(0, top_slot-1, F->linkNode(down));
                    F->createReducedNode(nkp, v, p);
                    assert(v.isVoid());

                    nk->setFull(i, p);
                    continue;
                }
                if (i == top_slot-1) {
                    unpacked_node* nkp = unpacked_node::newWritable(F,
                            -k, 1, SPARSE_ONLY);
                    nkp->setSparse(0, top_slot, F->linkNode(down));
                    F->createReducedNode(nkp, v, p);
                    assert(v.isVoid());

                    nk->setFull(i, p);
                    continue;
                }
                // no change
                if (F->isIdentityReduced()) {
                    nk->setFull(i, F->linkNode(down));
                } else {
                    unpacked_node* nkp = unpacked_node::newWritable(F,
                            -k, 1, SPARSE_ONLY);
                    nkp->setSparse(0, i, F->linkNode(down));
                    F->createReducedNode(nkp, v, p);
                    assert(v.isVoid());

                    nk->setFull(i, p);
                }
            } // for i

            F->createReducedNode(nk, v, p);
            assert(v.isVoid());
            F->unlinkNode(down);
            down = p;
        } // for k

        R.set(down);
        return true;
    }
}

// ******************************************************************
// *                       swapmodel  methods                       *
// ******************************************************************

swapmodel::swapmodel()
{
    N = -1;
    for_each_slot_which_number = true;
}

// ******************************************************************

bool swapmodel::process_arg(unsigned &i, const char** argv)
{
    if (rsgen_model::process_arg(i, argv)) return true;

    if (0==strcmp("-array", argv[i])) {
        ++i;
        for_each_slot_which_number = true;
        return true;
    }
    if (0==strcmp("-number", argv[i])) {
        ++i;
        for_each_slot_which_number = false;
        return true;
    }

    if ('-' == argv[i][0])  return false;

    // assume this is the size
    N = atoi(argv[i]);
    ++i;
    return true;
}

// ******************************************************************

void swapmodel::show_usage(std::ostream &out, const char* arg0)
{
    show_synopsis(out, arg0, "[options] N");

    show_model_info(out, "An array of N unique integers, with operations to swap adjacent integers. Initially, the integers are in order.");

    show_section(out, "SIZE PARAMETERS");
    show_switch(out, "N", "number of array slots and unique numbers");

    show_section(out, "ENCODING OPTIONS");
    show_switch(out, "-array", "System state is, for each array slot, which number is there. Slot #0 is the bottom DD variable, and slot #N-1 is the top DD variable.");
    show_switch(out, "-number", "System state is, for each number, which array slot contains it. Number 0 is the bottom DD variable, and number N-1 is the top DD variable.");

    rsgen_model::show_usage(out);
}

// ******************************************************************

bool swapmodel::build_model(MEDDLY::output &out)
{
    if (N<0) {
        std::cerr << "Missing or bad array size.\n";
        return false;
    }
    out << "Array swaps model, array has " << N << " slots\n";
    if (for_each_slot_which_number) {
        out << "Encoding using 'for each array slot, which number is there'\n";
    } else {
        out << "Encoding using 'for each number, which array slot is it in'\n";
    }
    out.flush();

    //
    // Initialize domain
    //
    int* sizes = new int[N];
    for (int i=N-1; i>=0; i--) sizes[i] = N;
    set_domain(MEDDLY::domain::createBottomUp(sizes, N));
    delete[] sizes;

    //
    // Add events
    //

#ifdef CHECK_SWAP_EVENT
    using namespace MEDDLY;

    rsgen_model::event* e = nullptr;
    if (for_each_slot_which_number) {
        e = new swap_adjacent<true>(N, N-1);
    } else {
        e = new swap_adjacent<false>(N, N-1);
    }

    forest* F = forest::create(get_domain(), RELATION, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL);

    dd_edge R(F);

    e->build_nsf(R);

    out << "Event edge is " << R << "\n";
    F->dump(out, ~0);

    return false;
#else
    if (for_each_slot_which_number) {
        for (int k=1; k<N; k++) {
            add_event(new swap_adjacent<true>(N, k));
        }
    } else {
        for (int k=1; k<N; k++) {
            add_event(new swap_adjacent<false>(N, k));
        }
    }
#endif

    return true;
}

// ******************************************************************

void swapmodel::build_initial(MEDDLY::minterm &init)
{
    for (int k=0; k<N; k++) {
        init.setVar(k+1, k);
    }
}

// ******************************************************************

bool swapmodel::known_cardinality(double &theory) const
{
    theory = 1.0;
    for (int i=2; i<=N; i++) {
        theory *= i;
    }
    return true;
}

// ******************************************************************

#ifdef HAVE_LIBGMP

bool swapmodel::known_cardinality(mpz_t theory) const
{
    mpz_set_ui(theory, 1);
    for (int i=2; i<=N; i++) {
        mpz_mul_ui(theory, theory, i);
    }
    return true;
}

#endif

// **********************************************************************
// main
// **********************************************************************

int main(int argc, const char** argv)
{
    swapmodel S;

    unsigned nArgs = (unsigned)argc;
    unsigned i=1;

    while (i < nArgs) {
        if (!S.process_arg(i, argv)) {
            S.show_usage(std::cerr, argv[0]);
            return 1;
        }
    }
    if (!S.has_dimensions()) {
        S.show_usage(std::cerr, argv[0]);
        return 1;
    }

    try {
        MEDDLY::initialize();
        MEDDLY::ostream_output out(std::cout);

        if (!S.build_model(out)) {
            S.show_usage(std::cerr, argv[0]);
            return 1;
        }
        S.generate(out, true, true);

        MEDDLY::cleanup();
        return 0;
    }
    catch (MEDDLY::error e) {
        std::cerr   << "\nCaught meddly error '" << e.getName()
                    << "'\n    thrown in " << e.getFile()
                    << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        std::cerr   << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    std::cerr << "\nSome other error?\n";
    return 4;

}


