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
  The interface is similar to STL iterators, except it is slightly
  faster to use the boolean operator than to compare with ```end()```,
  i.e., to use a loop of the form
  ```c++
  for (dd_edge::iterator I = e.begin(); I; ++I)
  ```
  You can pass a mask as a parameter to ```begin()```;
  this should be a pointer to a (const) ```minterm``` with values set to
    ```DONT_CARE``` for variables that are free to take any value,
    ```DONT_CHANGE``` for relation "to" variables that must equal
    their "from" counterpart, or a regular non-negative integer value
    to fix a variable.


* Old ```forest::getElement()``` replaced by ```dd_edge::getElement()```.

* Old ```forest::createEdge()``` for constant functions
  replaced by ```forest::createConstant()```.

### Implementation

