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

### Implementation

