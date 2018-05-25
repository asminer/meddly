---
title: Features
---

## Legend

| Symbol | Meaning |
| :----: | ------- |
| <span style="color:blue"> Y </span> | Implemented and stable |
| <span style="color:green"> y </span> | Implemented and unstable (untested, or known to have bugs, or likely to change) |
| <span style="color:orange"> n </span> | Under development |
| <span style="color:brown"> N </span> | Planned |
|   | No plans |

## Forest types

| Abbrv. | Status | Range | Value mechanism | Set/relation |
| :----- | :----: | :---- | :-------------- | :----------- |
| MDD | <span style="color:blue"> Y </span> | boolean | terminals | set |
| MxD | <span style="color:blue"> Y </span> | boolean | terminals | relation |
| MTMDD | <span style="color:blue"> Y </span> | integer | terminals | set |
| MTMxD | <span style="color:blue"> Y </span> | integer | terminals | relation |
| MTMDD | <span style="color:blue"> Y </span> | real | terminals | set |
| MTMxD | <span style="color:blue"> Y </span> | real | terminals | relation |
| MTMDD | <span style="color:brown"> N </span> | user-defined | terminals | set |
| MTMxD | <span style="color:brown"> N </span> | user-defined | terminals | relation |
| EV+MDD | <span style="color:green"> y </span> | integer | sum of edge values | set |
| EV+MxD | <span style="color:green"> y </span> | integer | sum of edge values | relation |
| EV*MDD | <span style="color:green"> y </span> | real | product of edge values | set |
| EV*MxD | <span style="color:green"> y </span> | real | product of edge values | relation |

## Operations

| Operation | MDD | MxD | MTMDD | MTMxD | EV+MDD | EV+MxD | EV*MDD | EV*MxD |
| --------- |:---:|:---:|:-----:|:-----:|:------:|:------:|:------:|:------:|
| Node count | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | 
| Edge count | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | 
| Create edge | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:green"> y </span> | <span style="color:green"> y </span> | <span style="color:green"> y </span> | <span style="color:green"> y </span> | 
| Iterators | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:green"> y </span> | <span style="color:green"> y </span> | <span style="color:green"> y </span> | <span style="color:green"> y </span> | 
| Cardinality | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:green"> y </span> | <span style="color:green"> y </span> | <span style="color:green"> y </span> | <span style="color:green"> y </span> | 
| Complement  | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> |
| Or / Set Union  | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> |
| And / Set Intersection  | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> |
| Set Difference  | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> |
| Pre-Image  | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Post-Image  | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Forward reachability | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Forward saturation | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Backward reachability | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Backward saturation | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Element-wise + | | | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Element-wise - | | | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Element-wise * | | | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Element-wise / | | | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Element-wise compare | | | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
| Element-wise max/min | | | <span style="color:blue"> Y </span> | <span style="color:blue"> Y </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> | <span style="color:orange"> n </span> |
