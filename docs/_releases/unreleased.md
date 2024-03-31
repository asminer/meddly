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

* For numerical operations, the arguments are
    whatever was passed previously to build the
    specialized operation, but without the ```argument```
    class hierarchy.

* For saturation operations, the arguments are the initial states
    forest, the relation, and the result forest.

As a result, the entire ```opname``` class hierarchy is now
obsolete, and has been removed.


### Deprecated methods

Uncomment ```ALLOW_DEPRECATED_0_17_5``` in ```defines.h```
to use these deprecated methods:

* ```destroyOperation```: use ```operation::destroy``` instead

* ```getOperation```: the builtin operation can be used as a function
    in the same way. For example, instead of ```getOperation(UNION, a, b, c)```,
    use ```UNION(a, b, c)```.

### Renamed classes

* ```satpregen_opname::pregen_relation``` is now ```pregen_relation```

* ```satotf_opname::subevent``` is now ```otf_subevent```
* ```satotf_opname::event``` is now ```otf_event```
* ```satotf_opname::relation``` is now ```otf_relation```

* ```satimpl_opname::implicit_relation``` is now ```implicit_relation```

### Implementation

