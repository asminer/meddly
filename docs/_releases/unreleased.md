---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Interface Changes

Builtin operations are now functions with arguments.

* For unary operations, the arguments are the forest to
    pull the input argument from, and either the forest
    for the output result, or the type of output result.

* For binary operations, the arguments are the forests
    for the input arguments and output result.

* For specialized operations, the arguments are
    whatever was passed previously to build the
    specialized operation, but without the ```argument```
    class hierarchy.

As a result, the entire ```opname``` class hierarchy is now
obsolete, and has been removed.


### Deprecated methods

Uncomment ```ALLOW_DEPRECATED_0_17_5``` in ```defines.h```
to use these deprecated methods:

* ```destroyOperation```: use ```operation::destroy``` instead

* ```getOperation```: the builtin operation can be used as a function
    in the same way. For example, instead of ```getOperation(UNION, a, b, c)```,
    use ```UNION(a, b, c)```.

### Implementation

