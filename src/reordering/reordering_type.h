// $Id$

#ifndef REORDERING_TYPE_H
#define REORDERING_TYPE_H

enum class reordering_type {
  // Always choose the lowest swappable inversion
  LOWEST_INVERSION,
  // Always choose the highest swappable inversion
  HIGHEST_INVERSION,
  // Sink down the "heaviest" variables
  SINK_DOWN,
  // Bubble up the "lightest" variables
  BUBBLE_UP,
  // Always choose the swappabel inversion where the variable on the top
  // has the fewest associated nodes
  LOWEST_COST,
  // Always choose the swappable inversion that will result in
  // the lowest memory consumption
  LOWEST_MEMORY,
  // Choose the swappable inversion randomly
  RANDOM,
  // Always choose the swappable inversion with the lowest average reference count
  LARC
};

#endif
