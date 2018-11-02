
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

namespace MEDDLY {
  class storage_initializer;
};

class MEDDLY::storage_initializer : public initializer_list {

    node_storage_style* simple;
    node_storage_style* pattern;
    node_storage_style* best;

    // TBD - build a new compact style

  public:
    storage_initializer(initializer_list *p);
  protected:
    virtual void setup();
    virtual void cleanup();
};

