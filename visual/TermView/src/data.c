
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

#include "data.h"

#include <stdlib.h>

void initialize(forest_t *f)
{
  if (0==f) return;
  f->fid = 0;
  f->name = 0;
  f->left = 0;
  f->right = 0;
  f->counts_raw = 0;
  f->counts = 0;
}

void destroy(forest_t *f)
{
  if (0==f) return;

  free(f->name);
  free(f->counts_raw);

  /* For sanity */
  f->name = 0;
  f->counts_raw = 0;
  f->counts = 0;
}


update_t* free_list = 0;

update_t* new_update(int fid, int level, int delta, update_t* next)
{
  update_t* node = 0;
  if (free_list) {
    node = free_list;
    free_list = free_list->next;
  } else {
    node = (update_t*) malloc(sizeof(update_t));
  }
  node->fid = fid;
  node->level = level;
  node->delta = delta;
  node->next = next;
  return node;
}

void kill_update(update_t* node)
{
  while (node) {
    update_t* nxt = node->next;
    node->next = free_list;
    free_list = node;
    node = nxt;
  }
}


