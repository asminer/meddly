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

* Added a ```relation forest``` class;
    the idea is to move implicit relation nodes here.
    Added support for binary operations between a (set) forest
    and a relation forest, with a result in a (set) forest;
    this will be used for relational product and saturation.

### Implementation

* Forest I/O methods are now centralized and non-virtual;
    all implementation removed from derived classes.

