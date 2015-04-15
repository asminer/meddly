
/* $Id$ */

/*
    termview: terminal viewing utility for Meddly trace data.
    Copyright (C) 2015, Iowa State University Research Foundation, Inc.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along 
    with this software.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DATA_H
#define DATA_H

/*
  Struct for forest information
*/
typedef struct {
  int fid;
  char* name;
  int left;
  int right;
  int* counts_raw;
  int* counts;
} forest_t;

/*
  Initialize a forest_t struct to zeroes.
*/
void initialize(forest_t *f);

/*
  Free the memory used by a forest_t struct.
  Does not free the pointer to the forest_t itself.
*/
void destroy(forest_t *f);

#endif
