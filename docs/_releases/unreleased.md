---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Interface Changes

In class ```dd_edge```:

* Added method ```showRelation()``` for better display of relations.

* Added (limited) ability to validate an MDD against a file.
  Essentially this is for debugging, in case there is ever
  an issue with a written MDD not matching a read MDD.


### Implementation

* Better node marking and displaying as needed for showRelation.

* Multi-terminal distance function change: unreachable states
  now have distance 0. Initial states now have distance 1.
