#!/bin/bash
#
#
# Massive script to run tests
#

EXDIR=examples

if [ ! -d $EXDIR ]; then
    EXDIR=../examples
fi
if [ ! -d $EXDIR ]; then
    echo "Can't find examples directory; run this script in"
    echo "the root directory or within developers/."
    exit 1
fi

print_status()
{
  if [ "x$2" != "x" ]; then
    printf "%20s %11s\n" $1 "$2"
  fi
}

summarize()
{
  echo "==============================================================="
  print_status "nqueens" "$s_nqueens"
  print_status "queen_cover" "$s_qcover"
  print_status "phils-bfs" "$s_philsbfs"
  print_status "kanban-bfs" "$s_kanbanbfs"
  print_status "slot-bfs" "$s_slotbfs"
  print_status "phils-dfs" "$s_philsdfs"
  print_status "kanban-dfs" "$s_kanbandfs"
  print_status "slot-dfs" "$s_slotdfs"
  print_status "kanban-exp" "$s_kanbanexp"
  print_status "slot-exp" "$s_slotexp"
  echo "==============================================================="
}

#
# Args: cmd
run_or_exit()
{
    if ! $@; then
        echo
        echo "Fail on: $@"
        echo
        exit 1
    fi
}

#
#  Super fast tests
#

for i in 1 2 3 4 5; do
  run_or_exit $EXDIR/nqueens $i
  s_nqueens="1..$i "
  summarize
done

for i in 1 2 3 4 5; do
  run_or_exit $EXDIR/queen_cover $i
  s_qcover="1..$i "
  summarize
done

for i in 2 4 6 8 10 15 20; do
  run_or_exit $EXDIR/dining_phils -n$i
  s_philsbfs="..$i "
  summarize
done

for i in 1 2 3 4 5 6; do
  run_or_exit $EXDIR/kanban $i -bfs
  s_kanbanbfs="..$i "
  summarize
done

for i in 2 3 4 5 6 7 8 9; do
  run_or_exit $EXDIR/slot $i -bfs
  s_slotbfs="2..$i "
  summarize
done

for i in 2 4 6 8 10 15 20; do
  run_or_exit $EXDIR/dining_phils -n$i -dfs
  s_philsdfs="..$i "
  summarize
done

for i in 1 2 3 4 5 6; do
  run_or_exit $EXDIR/kanban $i -dfs
  s_kanbandfs="..$i "
  summarize
done

for i in 2 4 6 8 10 15; do
  run_or_exit $EXDIR/slot $i -dfs
  s_slotdfs="..$i "
  summarize
done

for i in 1 2 3; do
  run_or_exit $EXDIR/kanban $i -exp
  s_kanbanexp="1..$i "
  summarize
done

for i in 2 3 4; do
  run_or_exit $EXDIR/slot $i -exp
  s_slotexp="2..$i "
  summarize
done

#
#  Medium tests
#

for i in 6 7 8 9 10 11 12; do
  run_or_exit $EXDIR/nqueens $i
  s_nqueens="1..$i+"
  summarize
  s_nqueens="1..$i "
done

for i in 6 7 8 9 10; do
  run_or_exit $EXDIR/queen_cover $i
  s_qcover="1..$i+"
  summarize
  s_qcover="1..$i "
done

for i in 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 200; do
  run_or_exit $EXDIR/dining_phils -n$i
  s_philsbfs="..$i+"
  summarize
  s_philsbfs="..$i "
done

for i in 7 8 9 10 15 20 25 30 35 40 45 50 55 60; do
  run_or_exit $EXDIR/kanban $i -bfs
  s_kanbanbfs="..$i+"
  summarize
  s_kanbanbfs="..$i "
done

for i in 10 11 12 13 14 15 16 17; do
  run_or_exit $EXDIR/slot $i -bfs
  s_slotbfs="2..$i+"
  summarize
  s_slotbfs="2..$i "
done

for i in 25 50 75 100 250 500 750 1000 1250 2500 5000; do
  run_or_exit $EXDIR/dining_phils -n$i -dfs
  s_philsdfs="..$i+"
  summarize
  s_philsdfs="..$i "
done

for i in 7 8 9 10 15 20 25 30 35 40 45 50 75 100 125; do
  run_or_exit $EXDIR/kanban $i -dfs
  s_kanbandfs="..$i+"
  summarize
  s_kanbandfs="..$i "
done

for i in 20 25 30 35 40 45; do
  run_or_exit $EXDIR/slot $i -dfs
  s_slotdfs="..$i+"
  summarize
  s_slotdfs="..$i "
done

for i in 4; do
  run_or_exit $EXDIR/kanban $i -exp
  s_kanbanexp="1..$i+"
  summarize
  s_kanbanexp="1..$i "
done

for i in 5; do
  run_or_exit $EXDIR/slot $i -exp
  s_slotexp="2..$i+"
  summarize
  s_slotexp="2..$i "
done

#
#  Long tests
#

for i in 13 14; do
  run_or_exit $EXDIR/nqueens $i
  s_nqueens="1..$i+"
  summarize
  s_nqueens="1..$i "
done

for i in 11 12; do
  run_or_exit $EXDIR/queen_cover $i
  s_qcover="1..$i+"
  summarize
  s_qcover="1..$i "
done

for i in 300 400 500; do
  run_or_exit $EXDIR/dining_phils -n$i
  s_philsbfs="..$i+"
  summarize
  s_philsbfs="..$i "
done

for i in 65 70 75; do
  run_or_exit $EXDIR/kanban $i -bfs
  s_kanbanbfs="..$i+"
  summarize
  s_kanbanbfs="..$i "
done

for i in 18 19 20; do
  run_or_exit $EXDIR/slot $i -bfs
  s_slotbfs="2..$i+"
  summarize
  s_slotbfs="2..$i "
done

for i in 6000 7000 8000; do
  run_or_exit $EXDIR/dining_phils -n$i -dfs
  s_philsdfs="..$i+"
  summarize
  s_philsdfs="..$i "
done

for i in 150 175 200; do
  run_or_exit $EXDIR/kanban $i -dfs
  s_kanbandfs="..$i+"
  summarize
  s_kanbandfs="..$i "
done

for i in 50 75 100; do
  run_or_exit $EXDIR/slot $i -dfs
  s_slotdfs="..$i+"
  summarize
  s_slotdfs="..$i "
done

for i in 5; do
  run_or_exit $EXDIR/kanban $i -exp
  s_kanbanexp="1..$i+"
  summarize
  s_kanbanexp="1..$i "
done

for i in 6; do
  run_or_exit $EXDIR/slot $i -exp
  s_slotexp="2..$i+"
  summarize
  s_slotexp="2..$i "
done


summarize
