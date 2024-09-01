---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Overview

Most changes this release are to clean up the compute table interface
and the operation interfaces,
to make it easier to implement operations.


### Interface Changes

* The ```unary_operation``` and ```binary_operation``` interfaces
  are changing; the old interface is still supported for now.
  The main differences are new virtual compute() methods that are "low level" for speed
  but also generic to minimize the number of overloaded compute() methods.
  See files ```oper_binary.h``` and ```oper_unary.h``` for details.

* The ```compute_table``` class now maintains a list of all compute
    tables, mainly for garbage collection purposes.

* ```operation::showAllComputeTables``` is now ```compute_table::showAll```

* An ```operation``` object may now be deleted directly,
  instead of using ```operation::destroy()```.

* Added new object, ```ct_itemtype```, for type-checking compute table entries.

* Existing object, ```ct_entry_type```, is now based on vectors of ```ct_itemtype``` objects.

* The registry of compute table entry types has been moved to
  class ```ct_entry_type```.

* The monolithic compute table has been moved from class ```operation```
  to class ```compute_table```.

* New method ```forest::isSingletonNode``` replaces old
    ```forest::getSingletonIndex``` and ```forest::getSingletonDown``` methods.

* New methods ```forest::makeRedundantsTo``` and ```forest::makeIdentitiesTo```
  should be used when converting between reduction types (typically
  in operations).

### Implementation

* New unified compute table implementation in file ```ct_styles.cc```,
  based on templates.

* Added possibility for huge (more than 4 billion entries) compute tables.
  This is a compute table setting, with a default of non-huge tables
  to save memory.

* Node reduction (in class ```forest```) is now a local operation;
  helper virtual methods (```normalize```, ```isIdentityEdge```, etc)
  have been removed.

* Boolean (MDD or MXD) Intersection, Union, Difference, and Complement operations
  have been rewritten using the new compute table and operation interfaces.

* Copy operations have been rewritten using the new compute table
  and operation interfaces.
