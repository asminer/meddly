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

#ifndef MEDDLY_PREPOSTIMAGE_H
#define MEDDLY_PREPOSTIMAGE_H

namespace MEDDLY {
    class binary_factory;
    binary_factory& PRE_IMAGE();
    binary_factory& POST_IMAGE();

    // TBD below here

    class forest;
    class binary_operation;

    /// The 'transitive closure postimage' operation builder.
    binary_operation* TC_POST_IMAGE(forest* a, forest* b, forest* c);
    void TC_POST_IMAGE_init();
    void TC_POST_IMAGE_done();

    /// The 'vector matrix multiplication' operation builder.
    binary_operation* VM_MULTIPLY(forest* a, forest* b, forest* c);
    void VM_MULTIPLY_init();
    void VM_MULTIPLY_done();

    /// The 'matrix vector multiplication' operation builder.
    binary_operation* MV_MULTIPLY(forest* a, forest* b, forest* c);
    void MV_MULTIPLY_init();
    void MV_MULTIPLY_done();
}

#endif

