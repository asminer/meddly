---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Changes to variables

 Variable names are now ```std::string``` objects instead
 of ```char *```.
For an un-named variable, use an empty string.

### Changes to operations
Builtin operations are now functions with arguments.

* For unary operations, the arguments are the forest to
    pull the input argument from, and either the forest
    for the output result, or the type of output result.

* For binary operations, the arguments are the forests
    for the input arguments and output result.

* For ternary operations (for example, constrained saturation),
    the arguments are the forests
    for the input arguments and output result.

* For numerical operations, the arguments are
    whatever was passed previously to build the
    specialized operation, but without the ```argument```
    class hierarchy.

* For saturation operations, the arguments are the initial states
    forest, the relation, and the result forest.

### Deprecated methods

Uncomment ```ALLOW_DEPRECATED_0_17_5``` in ```defines.h```
to use these deprecated methods:

* ```destroyOperation```: use ```operation::destroy``` instead

* ```getOperation```: the builtin operation can be used as a function
    in the same way. For example, instead of ```getOperation(UNION, a, b, c)```,
    use ```UNION(a, b, c)```,
    where ```a``` and ```b``` are the forests for the operands,
    and ```c``` is the forest for the result.

### Renamed classes

* ```satpregen_opname::pregen_relation``` is now ```pregen_relation```

* ```satotf_opname::subevent``` is now ```otf_subevent```
* ```satotf_opname::event``` is now ```otf_event```
* ```satotf_opname::relation``` is now ```otf_relation```

* ```satimpl_opname::implicit_relation``` is now ```implicit_relation```

* ```satimpl_opname::implicit_relation``` is now ```implicit_relation```

* ```sathyb_opname::subevent``` is now ```hybrid_subevent```
* ```sathyb_opname::event``` is now ```hybrid_event```
* ```sathyb_opname::hybrid_relation``` is now ```hybrid_relation```

### Removed classes

*  The entire ```opname``` class hierarchy is now deprecated,
    and has been removed.

*  Class ```specialized_operation``` is no longer needed,
    and has been removed.
