---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### New features

* Added application ```applyinfo``` in directory ```examples```,
    for displaying documentation for built-in operations.

* Added ```DIST_MIN``` binary operation, for minimum distance where
  negatives may be used for infinity.

* Added ```DIST_INC``` unary operation, for incrementing distances
  where negatives may be used for infinity.

* Added tests for pre/post image operations.

* Pre- and post-image operations work for both sets of states,
  and distance functions. Distance functions may be stored using
  EV+MDDs, with infinity for unreachable states, or using multi-terminal,
  integer range MDDs, with negative values for unreachable states.


### Interface Changes

* In the ```forest``` class, overloaded virtual methods
    ```createEdgeForVar``` have been deprecated in favor of
    a single ```createEdgeForVar``` method, where the terms
    are generic ```rangeval``` objects.
    Also there is now a single implementation in the ```forest``` class,
    instead of various methods in derived classes.

* Added a new ```rel_node``` class,
    for abstracting (read-only) relation access.
    Different forests may store relations differently, internally,
    and might not be based on decision diagram nodes.

* Added a ```unary_factory``` class, for built-in unary operations,
    with documentation.

* Added a ```binary_factory``` class, for built-in binary operations,
    with documentation.

* Edge values: some getters replaced by type conversion operators

* Forest class: older getDownPtr methods are now deprecated

### Implementation

* Forest I/O methods are now centralized and non-virtual;
    all implementation removed from derived classes.

