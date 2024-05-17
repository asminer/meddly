---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Overview

Most changes this release are to clean up the compute table interface,
to make it easier to implement operations.
TBD...

### Interface Changes

* The ```compute_table``` class now maintains a list of all compute
    tables, mainly for garbage collection purposes.

* ```operation::showAllComputeTables``` is now ```compute_table::showAll```


* Added new object, ```ct_itemtype```, for type-checking compute table entries.

* Existing object, ```ct_entry_type```, is now based on vectors of ```ct_itemtype``` objects.

* The registry of compute table entry types has been moved to
  class ```ct_entry_type```.

* The monolithic compute table has been moved from class ```operation```
  to class ```compute_table```.

### Implementation

* New unified compute table implementation in file ```ct_styles.cc```,
  based on templates.

* Added possibility for huge (more than 4 billion entries) compute tables.

