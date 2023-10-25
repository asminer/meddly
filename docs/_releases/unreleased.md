---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### New objects

* ```edge_value```: object for generic edge values.

* ```terminal```: object for terminal nodes.

### Interface changes

Mostly I/O related, but overall the interface will gradually replace
overloaded methods for specific edge value types, with a single method
using the generic edge value object.
Similar for terminal nodes.


