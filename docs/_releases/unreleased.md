---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

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

* Added a ```unary_factory``` class, for built-in unary operations.

### Implementation

* Forest I/O methods are now centralized and non-virtual;
    all implementation removed from derived classes.

