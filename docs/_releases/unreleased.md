---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### New features

* Users can build an arbitrary unary operation by writing a function
    of the form ```void F(const rangeval &in, rangeval &out)```,
    creating an instance of object ```user_unary_factory```,
    and then using ```apply``` on that object.
    See the test file ```tests/ops_user_un.cc```.

### Interface Changes

* Built-in operations ```REACHABLE_STATES_BFS```
    and ```REVERSE_REACHABLE_BFS```
    are deprecated;
    use instead ```REACHABLE_TRAD_NOFS(true)```
    and ```REACHABLE_TRAD_NOFS(false)```, respectively.

* To read a domain from a file,
    method ```domain::read()``` has been replaced by
    an overloaded static method ```domain::create()```.

* Added "options" to operations.
    This is the mechanism to adjust algorithms, for example,
    different versions of saturation. Options are used only for
    adjustments that do not affect the final answer.

### Implementation

* The ```COPY``` operation will use relation nodes in some cases.
  This will allow conversions from implicit representations
  once those implementations are migrated over.
