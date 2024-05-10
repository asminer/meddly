---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Interface Changes

* Added new object, ```ct_itemtype```, for type-checking compute table entries.

* Existing object, ```ct_entry_type```, is now based on vectors of ```ct_itemtype``` objects.


### Implementation

* New unified compute table implementation in file ```ct_styles.cc```,
  based on templates.

* Added possibility for huge (more than 4 billion entries) compute tables.

