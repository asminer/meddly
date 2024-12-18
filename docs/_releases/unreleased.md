---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Interface Changes

* Added classes ```minterm``` and ```minterm_coll``` to deal
    with variable assignments and truth tables over sets and relations.

* Can create functions from single minterms or minterm collections
    using methods ```minterm::buildFunction()``` and
    ```minterm_coll::buildFunction```.
    This replaces several ```forest::createEdge()``` methods,
    which are now deprecated. (TBD)

* Function evaluation is now done through ```dd_edge::evaluate```.
    TBD: deprecate the old function evaluate

* ```dd_edge::set()``` parameters are reversed, for consistency
    (edge value before node).

* Iterating through a function is now done using
    an instance of ```dd_edge::iterator```, created with
    ```dd_edge::begin()```.

* Old ```forest::getElement()``` replaced by ```dd_edge::getElement()```.

### Implementation

