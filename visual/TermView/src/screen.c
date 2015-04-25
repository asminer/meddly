
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

#include "screen.h"
#include "system.h"

int open_screen(screen_t *S)
{
  if (0==initscr()) return 0;
  // check screen size
  int row, col;
  getmaxyx(stdscr, row, col);
  if (S) {
    S->rows = row;
    S->cols = col;
    S->forwidth = col;
  }
  if (col < 40) {
    endwin();
    fprintf(stderr, "\nTerminal width of at least 40 columns is required\n");
    return 0;
  }
  if (row < 10) {
    endwin();
    fprintf(stderr, "\nTerminal height of at lest 10 rows is required\n");
    return 0;
  }
  noecho();
  cbreak();
  clear();
  return 1;
}

void close_screen()
{
  endwin();
}

void update_keys(screen_t *S)
{
  if (0==S) return;
  move(S->rows-3, -2 + S->cols / 4);
  addch('P' | A_BOLD);
  if (S->playing) {
    addch('a');
    addch('u');
    addch('s');
    addch('e');
  } else {
    addch('l');
    addch('a');
    addch('y');
    addch(' ');
  }
  move(S->rows-3, -2+ 3*S->cols / 4);
  addch('Q' | A_BOLD);
  addch('u');
  addch('i');
  addch('t');
  move(S->rows-1, S->cols-1);
  refresh();
}

inline int level2row(screen_t *S, int level)
{
  if (0==level) return -1;

  if (S->has_primed) {
    if (level<0) {
      level *= -2;
      level--;
    } else {
      level *= 2;
    }
  } 

  return S->rows - 4 - level;
}

inline void drawValue(screen_t *S, int fid, int level, int value, int blanks)
{
  int row = level2row(S, level);
  int col = (fid-1) * S->forwidth + S->varwidth + 6;
  if (blanks) {
    move(row, col);
    for (int spaces = S->forwidth - (S->varwidth+6); spaces; spaces--)
      addch(' ');
  }
  mvprintw(row, col, "%d", value);
}

void update_windows(screen_t *S)
{
  if (0==S) return;

  /*
    Draw all the lines
  */
  int nb;
  for (int r=0; r<S->rows-4; r++) {
    move(r, 0);
    nb = S->nf;
    for (int c=0; c<S->cols; c++) {
      if (nb && (0== c % S->forwidth)) {
        addch(ACS_VLINE);
        nb--;
      } else {
        addch(' ');
      }
    }
  }
  move(S->rows-4, 0);
  if (S->nf > 0) {
    addch(ACS_LTEE);
  } else {
    addch(ACS_ULCORNER);
  }
  nb = S->nf;
  for (int c=1; c<S->cols-1; c++) {
    if (nb && (0== c % S->forwidth)) {
      addch(ACS_BTEE);
      nb--;
    } else {
      addch(ACS_HLINE);
    }
  }
  addch(ACS_URCORNER);
  move(S->rows-3, 0);
  addch(ACS_VLINE);
  for (int c=1; c<S->cols-1; c++) addch(' ');
  addch(ACS_VLINE);
  move(S->rows-2, 0);
  addch(ACS_LLCORNER);
  for (int c=1; c<S->cols-1; c++) addch(ACS_HLINE);
  addch(ACS_LRCORNER);

  /*
    Draw the variables
  */
  for (int f=0; f<S->maxf; f++) {
    forest_t* F = S->F[f];
    if (0==F) continue;
    int c = (f-1) * S->forwidth + 2;
    for (int k=F->left; k<=F->right; k++) {
      int r = level2row(S, k);
      if (r<0) continue;
      move(r, c);
      if (F->left<0) {
        addch('x');
        if (k<0) {
          printw("'%0*d: ", S->varwidth, -k);
        } else {
          printw(" %0*d: ", S->varwidth, k);
        }
      } else {
        addch(' ');
        addch('x');
        printw("%0*d: ", S->varwidth, k);
      }
      drawValue(S, F->fid, k, F->counts[k], 0);
    }
    
  }
  
  /*
    Keystroke area
  */
  update_keys(S);
}

void update_p(const char* str)
{
  int row, col;
  getmaxyx(stdscr, row, col);
  move(row-1, 0);
  int i;
  for (i=0; i<col; i++) {
    if (0==str[i]) break;
    addch(str[i]);
  }
  for (; i<col; i++) {
    addch(' ');
  }
  refresh();
}

int keystroke(char block)
{
  if (block) {
    timeout(-1);
  } else {
    timeout(0);
  }
  return getch();
}

void update_a(screen_t* S, update_t* list)
{
  if (0==S) return;
  while (list) {
    forest_t* F = 0;
    if (list->fid > 0 && list->fid < S->maxf) {
      F = S->F[list->fid];
    }
    if (F) if (list->level >= F->left && list->level <= F->right) {
      F->counts[list->level] += list->delta;
      drawValue(S, list->fid, list->level, F->counts[list->level],
        (list->delta < 0) ? 1 : 0);
    }

    list = list->next;
  }
  move(S->rows-1, S->cols-1);
  refresh();
}