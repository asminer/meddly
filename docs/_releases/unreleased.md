---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Interface Changes

* Added classes ```minterm``` and ```minterm_coll``` to deal
    with variable assignments and truth tables over sets and relations.

* Added operators for MDDs and collections of minterms.
    TBD: deprecate old function::createEdge that takes
    in a list of minterms.

* Function evaluation is now done through ```dd_edge::evaluate```.
    TBD: deprecate the old function evaluate

* ```dd_edge::set()``` parameters are reversed, for consistency
    (edge value before node).

### Implementation

