---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Interface Changes

Most changes are simplifications and reorganization of class ```dd_edge```.
Some method changes:

* Use method ```attach()``` in place of ```setForest()```.

* For method ```writePicture```, using an extension of ```dot```
  will produce a dot file, without running any of the Graphviz utilities.

* Method ```show()``` does not take the verbosity level parameter
  any more.
  To display the graph, use new method ```showGraph'''.

* Method ```getCardinality()``` has been removed
    (use apply(CARDINALITY...) instead).

* Operators are still supported, but the interface and implementation
  have been moved to operators.h/.cc.

