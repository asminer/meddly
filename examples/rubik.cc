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
        /// Rotation event base class
        class rotate : public rsgen_model::event {
            protected:
                const std::vector <int> &order;
                char face;  // U, D, F, B, L, R
                int turns;  // 1: clockwise, 2: flip, 3: counter clockwise
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
                        @param  rotdir      For each variable k, rotdir[k]
                                            gives the direction of rotation as
                                            follows:
                                               'i': no change
                                               'x': swap y,z
                                               'y': swap x,z
                                               'z': swap x,y
                                            anything else is no constraint.
                        @param  r           Edge to store the relation.
                */
                void xp_rotate_y(int xp_lvl, int y_lvl,
                        const std::vector <char> &rotdir, dd_edge &r) const;

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
        const std::vector <char> &rotdir, dd_edge &r) const
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
                if ('x' == rotdir[k]) {
                    m.setVars(k, DONT_CARE, (0==y) ? 0 : 3-y);  // swap 2,3
                    continue;
                }
                if ('y' == rotdir[k]) {
                    m.setVars(k, DONT_CARE, (1==y) ? 1 : 2-y);  // swap 0,2
                    continue;
                }
                if ('z' == rotdir[k]) {
                    m.setVars(k, DONT_CARE, (2==y) ? 2 : 1-y);  // swap 0,1
                    continue;
                }
                assert(false);
                continue;
            } // u_x_lvl == k

            if (u_y_lvl == k) {
                m.setVars(k, y, DONT_CARE);
                continue;
            }

            if ('i' == rotdir[k]) {
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
    return false;
}

// *************************************************************
bool rubik::pos_rotate::build_nsf(MEDDLY::dd_edge &R)
{
    return false;
}

// ******************************************************************
// *                              main                              *
// ******************************************************************

int main()
{
    return 0;
}

