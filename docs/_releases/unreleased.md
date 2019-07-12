---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Simple Interface Changes

* Under forest policies, you can now decide (at runtime) to use
  reference counts or mark and sweep for garbage collection.
  The current default is to use reference counts.
  Set ```useReferenceCounts``` to false to instead use mark and sweep.

### Expert Interface Changes

* More compact (and dynamic) node header implementation.

### Implementation


