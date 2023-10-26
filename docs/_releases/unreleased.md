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

Source files ```encoders.h``` and ```encoders.cc``` were removed,
as this functionality is now captured by the new ```edge_value```
and ```terminal``` objects.

### Implementation

Note that the exchange format for reading/writing DDs has changed.
Files written with earlier versions of the library will not be readable.

