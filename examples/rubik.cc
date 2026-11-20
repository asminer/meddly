/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

/*
    Another Rubik's cube model.

    We assume the cube is oriented with its center at
    the origin of an x,y,z coordinate system with
        Up face having +z normal vector (usually white)
        Down face having -z normal vector (usually yellow)
        Front face having +x normal vector (usually blue)
        Back face having -x normal vector (usually green)
        Right face having +y normal vector (usually orange)
        Left face having -y normal vector (usually red)
    The entire cube is never rotated or moved.

    We name positions and pieces the same way,
    with color pairs (for edge pieces) or color triples (for corners).
    For the triples, the order is given by the faces read off
    in x,y,z normal vector order, while the cube is oriented so that
    the piece is in the +x,+y,+z orthant.

    There are two parts of the state. One is, which pieces are where.
    The second is, the orientation of each piece.

    For which pieces are where, we consider two different encodings.
        * position based: for each position, which piece is there
        * piece based   : for each piece, where is it

    Orientations are x, y, z,
        indicating which vector is normal to the first color
        in the piece name.


    Corner pieces and solved orientations:
        ULF, z
        UFR, z
        URB, z
        UBL, z
        DLB, z
        DBR, z
        DRF, z
        DFL, z

    Edge pieces and solved orientations:
        UL, z
        UF, z
        UR, z
        UB, z
        RB, y
        RF, y
        LB, y
        LF, y
        DL, z
        DF, z
        DR, z
        DB, z

    Moves are named by the face name as looking directly at the face.
    Names are
        U, for clockwise rotation of Up face
        U', for counter clockwise rotation of Up face
        U2, for Up flip (2 clockwise or counter rotations)
    etc.

    Up rotations:
        cycle pieces appropriately in positions ULF, UFR, URB, UBL
        cycle pieces appropriately in positions UL, UF, UR, UB
        all orientations:
            x -> y
            y -> x
            z -> z
        (no change for U2)

    Right rotations:
        cycle pieces appropriately in positions UFR, URB, DBR, DRF
        cycle pieces appropriately in positions UR, RB, RF, DR
        all orientations:
            x -> z
            y -> y
            z -> x
        (no change for R2)

    Front rotations:
        cycle pieces appropriately in positions ULF, UFR, DRF, DFL
        cycle pieces appropriately in positions UF, RF, LF, DF
        all orientations:
            x -> x
            y -> z
            z -> y
        (no change for F2)

 */

#include "meddly/meddly.h"
#include "utils/rsgen_model.h"

using namespace MEDDLY;

/*
    Corner piece names
*/
/*
const unsigned ULF_piece = 0;
const unsigned UFR_piece = 1;
const unsigned URB_piece = 2;
const unsigned UBL_piece = 3;
const unsigned DLB_piece = 4;
const unsigned DBR_piece = 5;
const unsigned DRF_piece = 6;
const unsigned DFL_piece = 7;
*/

/*
    Edge piece names
*/
/*
const unsigned UL_piece = 0;
const unsigned UF_piece = 1;
const unsigned UR_piece = 2;
const unsigned UB_piece = 3;
const unsigned RB_piece = 4;
const unsigned RF_piece = 5;
const unsigned LB_piece = 6;
const unsigned LF_piece = 7;
const unsigned DL_piece = 8;
const unsigned DF_piece = 9;
const unsigned DR_piece = 10;
const unsigned DB_piece = 11;
*/

/*
    Variable names.
    The variable ordering array should map these to levels.
    A level of 0 means that variable is unused in the model
    (e.g., 2x2 cubes, or disabling features).
*/
const unsigned var_ULF_piece    = 0;
const unsigned var_UFR_piece    = 1;
const unsigned var_URB_piece    = 2;
const unsigned var_UBL_piece    = 3;
const unsigned var_DLB_piece    = 4;
const unsigned var_DBR_piece    = 5;
const unsigned var_DRF_piece    = 6;
const unsigned var_DFL_piece    = 7;
const unsigned var_UL_piece     = 8;
const unsigned var_UF_piece     = 9;
const unsigned var_UR_piece     = 10;
const unsigned var_UB_piece     = 11;
const unsigned var_RB_piece     = 12;
const unsigned var_RF_piece     = 13;
const unsigned var_LB_piece     = 14;
const unsigned var_LF_piece     = 15;
const unsigned var_DL_piece     = 16;
const unsigned var_DF_piece     = 17;
const unsigned var_DR_piece     = 18;
const unsigned var_DB_piece     = 19;
const unsigned var_ULF_orient   = 20;
const unsigned var_UFR_orient   = 21;
const unsigned var_URB_orient   = 22;
const unsigned var_UBL_orient   = 23;
const unsigned var_DLB_orient   = 24;
const unsigned var_DBR_orient   = 25;
const unsigned var_DRF_orient   = 26;
const unsigned var_DFL_orient   = 27;
const unsigned var_UL_orient    = 28;
const unsigned var_UF_orient    = 29;
const unsigned var_UR_orient    = 30;
const unsigned var_UB_orient    = 31;
const unsigned var_RB_orient    = 32;
const unsigned var_RF_orient    = 33;
const unsigned var_LB_orient    = 34;
const unsigned var_LF_orient    = 35;
const unsigned var_DL_orient    = 36;
const unsigned var_DF_orient    = 37;
const unsigned var_DR_orient    = 38;
const unsigned var_DB_orient    = 39;


class rubik : public rsgen_model {
    public:
        /// Rotation event base class
        class rotate : public rsgen_model::event {
            protected:
                const std::vector <int> &order;
                char face;  // U, D, F, B, L, R
                int turns;  // 1: clockwise, 2: flip, 3: counter clockwise
                char face_axis; // 'x', 'y', or 'z'
            public:
                rotate(const char* name, char face, int turns,
                        const std::vector <int> &ord);

            protected:
                /**
                    Build (piece of) a relation with x' = y.
                    Intended for moving pieces in positional encoding.
                        @param  xp_lvl      Level for x'
                        @param  y_lvl       Level for y
                        @param  nochange    If nochange[k] then assume
                                            x'_k = x_k; otherwise no
                                            constraints on x'_k or x_k.
                        @param  r           Edge to store the relation.
                 */
                void xp_equals_y(int xp_lvl, int y_lvl,
                        const std::vector <bool> &nochange, dd_edge &r) const;

                /**
                    Build a relation with x' = rotate(y).
                    Intended for changing piece orientations.
                        @param  xp_lvl      Level for x'
                        @param  y_lvl       Level for y
                        @param  nochange    If nochange[k] then assume
                                            x'_k = x_k; otherwise no
                                            constraints on x'_k or x_k.
                        @param  r           Edge to store the relation.
                */
                void xp_rotate_y(int xp_lvl, int y_lvl,
                        const std::vector <bool> &nochange, dd_edge &r) const;


                inline int reorient(int v) const {
                    switch (face_axis) {
                        case 'x':
                            return (0==v) ? 0 : 3-v;    // swap 2,3

                        case 'y':
                            return (1==v) ? 1 : 2-v;    // swap 0,2

                        case 'z':
                            return (2==v) ? 2 : 1-v;    // swap 0,1

                        default:
                            assert(false);
                    }
                }

        };  // class rotate

        /// Rotations for encoding: for each piece, where is it
        class pos_rotate : public rotate {
                // 1 rotation will do A -> B -> C -> D -> A
                // 2 rotations will do A -> C -> A; B -> D -> B
                // 3 rotations will do A -> D -> C -> B -> A

                // Corner pieces (actual levels, or 0 if disabled)
                int cpA, cpB, cpC, cpD;
                // Corner orientations (actual levels, or 0 if disabled)
                int coA, coB, coC, coD;
                // Edge pieces (actual levels, or 0 if disabled)
                int epA, epB, epC, epD;
                // Edge orientations (actual levels, or 0 if disabled)
                int eoA, eoB, eoC, eoD;
            public:
                pos_rotate(const char* name, char face, int turns,
                        const std::vector <int> &ord);

                virtual bool fire(const MEDDLY::minterm &curr,
                        MEDDLY::minterm &);

                virtual bool build_nsf(MEDDLY::dd_edge &R);

            private:
                static void rotate_pieces(int A, int B, int C, int D,
                        minterm &x);
                void rotate_orientations(int A, int B, int C, int D,
                        minterm &x);

                static inline void rotate_levels(int &A, int &B, int &C, int &D)
                {
                    int t = D;
                    D = C;
                    C = B;
                    B = A;
                    A = t;
                }
        };


};

// ******************************************************************
// *                     rubik::rotate  methods                     *
// ******************************************************************

rubik::rotate::rotate(const char* name, char _face, int _turns,
        const std::vector <int> &ord) : event(name), order(ord)
{
    face = _face;
    turns = _turns;
    assert(turns>0);
    assert(turns<4);

    switch (face) {
        case 'F':
        case 'B':
                    face_axis = 'x';
                    break;

        case 'L':
        case 'R':
                    face_axis = 'y';
                    break;

        case 'U':
        case 'D':
                    face_axis = 'z';
                    break;

        default:
            face_axis = 0;
    }
}

// *************************************************************

void rubik::rotate::xp_equals_y(int xp_lvl, int y_lvl,
        const std::vector <bool> &nochange, dd_edge &r) const
{
    assert(xp_lvl < 0);
    assert(y_lvl > 0);

    const unsigned u_y_lvl = unsigned(y_lvl);
    const unsigned u_x_lvl = unsigned(-xp_lvl);

    forest* F = r.getForest();
    unsigned size = F->getVariableSize(y_lvl);

    minterm_coll mtc(size, F);

    for (unsigned y=0; y<size; y++) {
        minterm &m = mtc.unused();
        for (unsigned k=1; k<mtc.getNumVars(); k++) {
            if (nochange[k]) {
                // identity
                if (u_x_lvl == k) {
                    m.setVars(k, y, y);
                    continue;
                }
                if (u_y_lvl == k) {
                    m.setVars(k, y, y);
                    continue;
                }
                m.setVars(k, DONT_CARE, DONT_CHANGE);
            } else {
                // fully
                if (u_x_lvl == k) {
                    m.setVars(k, DONT_CARE, y);
                    continue;
                }
                if (u_y_lvl == k) {
                    m.setVars(k, y, DONT_CARE);
                    continue;
                }
                m.setVars(k, DONT_CARE, DONT_CARE);
            }
        } // for k
        m.setValue(true);
        mtc.pushUnused();
    }
    mtc.buildFunctionMax(false, r);
}

// *************************************************************

void rubik::rotate::xp_rotate_y(int xp_lvl, int y_lvl,
        const std::vector <bool> &nochange, dd_edge &r) const
{
    assert(xp_lvl < 0);
    assert(y_lvl > 0);

    const unsigned u_y_lvl = unsigned(y_lvl);
    const unsigned u_x_lvl = unsigned(-xp_lvl);

    forest* F = r.getForest();
    unsigned size = F->getVariableSize(y_lvl);
    assert(3==size);

    minterm_coll mtc(size, F);

    for (unsigned y=0; y<size; y++) {
        minterm &m = mtc.unused();
        for (unsigned k=1; k<mtc.getNumVars(); k++) {
            if (u_x_lvl == k) {
                m.setVars(k, DONT_CARE, reorient(y));
                continue;
            } // u_x_lvl == k

            if (u_y_lvl == k) {
                m.setVars(k, y, DONT_CARE);
                continue;
            }

            if (nochange[k]) {
                m.setVars(k, DONT_CARE, DONT_CHANGE);
            } else {
                m.setVars(k, DONT_CARE, DONT_CARE);
            }
        } // for k
        m.setValue(true);
        mtc.pushUnused();
    }
    mtc.buildFunctionMax(false, r);
}

// ******************************************************************
// *                   rubik::pos_rotate  methods                   *
// ******************************************************************

rubik::pos_rotate::pos_rotate(const char* name, char _face, int _turns,
        const std::vector <int> &ord) : rotate(name, _face, _turns, ord)
{
    // Set the pieces to be rotated
    switch (face) {
        case 'U':
            // corner pieces
            cpA = order[var_ULF_piece];
            cpB = order[var_UBL_piece];
            cpC = order[var_URB_piece];
            cpD = order[var_UFR_piece];
            // corner orientations
            coA = order[var_ULF_orient];
            coB = order[var_UBL_orient];
            coC = order[var_URB_orient];
            coD = order[var_UFR_orient];
            // edge pieces
            epA = order[var_UL_piece];
            epB = order[var_UB_piece];
            epC = order[var_UR_piece];
            epD = order[var_UF_piece];
            // edge orientations
            eoA = order[var_UL_orient];
            eoB = order[var_UB_orient];
            eoC = order[var_UR_orient];
            eoD = order[var_UF_orient];
            return;

        case 'D':
            // corner pieces
            cpA = order[var_DFL_piece];
            cpB = order[var_DRF_piece];
            cpC = order[var_DBR_piece];
            cpD = order[var_DLB_piece];
            // corner orientations
            coA = order[var_DFL_orient];
            coB = order[var_DRF_orient];
            coC = order[var_DBR_orient];
            coD = order[var_DLB_orient];
            // edge pieces
            epA = order[var_DL_piece];
            epB = order[var_DF_piece];
            epC = order[var_DR_piece];
            epD = order[var_DB_piece];
            // edge orientations
            eoA = order[var_DL_orient];
            eoB = order[var_DF_orient];
            eoC = order[var_DR_orient];
            eoD = order[var_DB_orient];
            return;

        case 'R':
            // corner pieces
            cpA = order[var_UFR_piece];
            cpB = order[var_URB_piece];
            cpC = order[var_DBR_piece];
            cpD = order[var_DRF_piece];
            // corner orientations
            coA = order[var_UFR_orient];
            coB = order[var_URB_orient];
            coC = order[var_DBR_orient];
            coD = order[var_DRF_orient];
            // edge pieces
            epA = order[var_UR_piece];
            epB = order[var_RB_piece];
            epC = order[var_DR_piece];
            epD = order[var_RF_piece];
            // edge orientations
            eoA = order[var_UL_orient];
            eoB = order[var_LB_orient];
            eoC = order[var_DL_orient];
            eoD = order[var_LF_orient];
            return;

        case 'L':
            // corner pieces
            cpA = order[var_ULF_piece];
            cpB = order[var_DFL_piece];
            cpC = order[var_DLB_piece];
            cpD = order[var_UBL_piece];
            // corner orientations
            coA = order[var_ULF_orient];
            coB = order[var_DFL_orient];
            coC = order[var_DLB_orient];
            coD = order[var_UBL_orient];
            // edge pieces
            epA = order[var_UL_piece];
            epB = order[var_LF_piece];
            epC = order[var_DL_piece];
            epD = order[var_LB_piece];
            // edge orientations
            eoA = order[var_UL_orient];
            eoB = order[var_LF_orient];
            eoC = order[var_DL_orient];
            eoD = order[var_LB_orient];
            return;

        case 'F':
            // corner pieces
            cpA = order[var_UFR_piece];
            cpB = order[var_DRF_piece];
            cpC = order[var_DFL_piece];
            cpD = order[var_ULF_piece];
            // corner orientations
            coA = order[var_UFR_orient];
            coB = order[var_DRF_orient];
            coC = order[var_DFL_orient];
            coD = order[var_ULF_orient];
            // edge pieces
            epA = order[var_UF_piece];
            epB = order[var_RF_piece];
            epC = order[var_DF_piece];
            epD = order[var_LF_piece];
            // edge orientations
            eoA = order[var_UF_orient];
            eoB = order[var_RF_orient];
            eoC = order[var_DF_orient];
            eoD = order[var_LF_orient];
            return;

        case 'B':
            // corner pieces
            cpA = order[var_URB_piece];
            cpB = order[var_UBL_piece];
            cpC = order[var_DLB_piece];
            cpD = order[var_DBR_piece];
            // corner orientations
            coA = order[var_URB_orient];
            coB = order[var_UBL_orient];
            coC = order[var_DLB_orient];
            coD = order[var_DBR_orient];
            // edge pieces
            epA = order[var_UB_piece];
            epB = order[var_LB_piece];
            epC = order[var_DB_piece];
            epD = order[var_RB_piece];
            // edge orientations
            eoA = order[var_UB_orient];
            eoB = order[var_LB_orient];
            eoC = order[var_DB_orient];
            eoD = order[var_RB_orient];
            return;

        default:
            assert(false);
    }
}

// *************************************************************

bool rubik::pos_rotate::fire(const MEDDLY::minterm &curr,
                        MEDDLY::minterm &next)
{
    next.setFrom(curr);
    for (int t=0; t<turns; t++) {
        rotate_pieces(cpA, cpB, cpC, cpD, next);
        rotate_pieces(epA, epB, epC, epD, next);
        rotate_orientations(coA, coB, coC, coD, next);
        rotate_orientations(eoA, eoB, eoC, eoD, next);
    }
    return true;
}

// *************************************************************

bool rubik::pos_rotate::build_nsf(MEDDLY::dd_edge &R)
{
    int ncpA = -cpA;
    int ncpB = -cpB;
    int ncpC = -cpC;
    int ncpD = -cpD;
    int ncoA = -coA;
    int ncoB = -coB;
    int ncoC = -coC;
    int ncoD = -coD;
    int nepA = -epA;
    int nepB = -epB;
    int nepC = -epC;
    int nepD = -epD;
    int neoA = -eoA;
    int neoB = -eoB;
    int neoC = -eoC;
    int neoD = -eoD;
    for (int t=0; t<turns; t++) {
        rotate_levels(ncpA, ncpB, ncpC, ncpD);
        rotate_levels(ncoA, ncoB, ncoC, ncoD);
        rotate_levels(nepA, nepB, nepC, nepD);
        rotate_levels(neoA, neoB, neoC, neoD);
    }

    //
    // build don't change set
    //
    std::vector <bool> identity(41);
    for (unsigned i=0; i<identity.size(); i++) {
        identity[i] = true;
    }
    identity[cpA] = false;
    identity[cpB] = false;
    identity[cpC] = false;
    identity[cpD] = false;
    identity[coA] = false;
    identity[coB] = false;
    identity[coC] = false;
    identity[coD] = false;
    identity[epA] = false;
    identity[epB] = false;
    identity[epC] = false;
    identity[epD] = false;
    identity[eoA] = false;
    identity[eoB] = false;
    identity[eoC] = false;
    identity[eoD] = false;

    //
    // Build portions of the relation, AND them together
    //
    dd_edge r(R);
    R.getForest()->createConstant(true, R);

    if (cpA && cpB && cpC && cpD)
    {
        xp_equals_y(ncpA, cpA, identity, r);
        R *= r;
        xp_equals_y(ncpB, cpB, identity, r);
        R *= r;
        xp_equals_y(ncpC, cpC, identity, r);
        R *= r;
        xp_equals_y(ncpD, cpD, identity, r);
        R *= r;
    }

    if (coA && coB && coC && coD)
    {
        xp_rotate_y(ncoA, coA, identity, r);
        R *= r;
        xp_rotate_y(ncoB, coB, identity, r);
        R *= r;
        xp_rotate_y(ncoC, coC, identity, r);
        R *= r;
        xp_rotate_y(ncoD, coD, identity, r);
        R *= r;
    }

    if (epA && epB && epC && epD)
    {
        xp_equals_y(nepA, epA, identity, r);
        R *= r;
        xp_equals_y(nepB, epB, identity, r);
        R *= r;
        xp_equals_y(nepC, epC, identity, r);
        R *= r;
        xp_equals_y(nepD, epD, identity, r);
        R *= r;
    }

    if (eoA && eoB && eoC && eoD)
    {
        xp_rotate_y(neoA, eoA, identity, r);
        R *= r;
        xp_rotate_y(neoB, eoB, identity, r);
        R *= r;
        xp_rotate_y(neoC, eoC, identity, r);
        R *= r;
        xp_rotate_y(neoD, eoD, identity, r);
        R *= r;
    }

    return true;
}

// *************************************************************

void rubik::pos_rotate::rotate_pieces(int A, int B, int C, int D, minterm &x)
{
    int tD = x.getVar(D);
    x.setVar(D, x.getVar(C));
    x.setVar(C, x.getVar(B));
    x.setVar(B, x.getVar(A));
    x.setVar(A, tD);
}

// *************************************************************

void rubik::pos_rotate::rotate_orientations(int A, int B, int C, int D,
        minterm &x)
{
    int tD = x.getVar(D);
    x.setVar(D, reorient(x.getVar(C)));
    x.setVar(C, reorient(x.getVar(B)));
    x.setVar(B, reorient(x.getVar(A)));
    x.setVar(A, reorient(tD));
}

// ******************************************************************
// *                              main                              *
// ******************************************************************

int main()
{
    try {

        std::vector <int> varorder(40);
        for (int i=0; i<40; i++) {
            varorder[i] = (i<8) ? i+1 : 0;
        }

        rubik::pos_rotate R("U", 'U', 3, varorder);

        MEDDLY::initialize();

        int sizes[40];
        for (int i=0; i<8; i++) sizes[i] = 8;
        for (int i=8; i<20; i++) sizes[i] = 12;
        for (int i=20; i<40; i++) sizes[i] = 3;
        domain* D = domain::createBottomUp(sizes, 40);

        forest* mxd = forest::create(D, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);

        dd_edge e(mxd);
        R.build_nsf(e);


        ostream_output out(std::cout);
        out << "Build NSF for event " << R.get_name() << ":\n";
        e.showRelation(out);


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


    return 0;
}

