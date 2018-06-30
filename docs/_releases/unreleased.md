---
title: New since last release
number: changes
layout: single
---

### Simple Interface Changes

* For forests with range of ```INTEGER```,
  function values are now ```long``` instead of ```int```.
  The interface has been updated appropriately
  (methods ```createEdge```, ```createEdgeForVar```, and ```evaluate```).

### Simple Interface Changes

 * Class ```compute_table::search_key```
   has been renamed ```compute_table::entry_key```, and is now concrete.

 * Class ```compute_table::search_result```
   has been renamed ```compute_table::entry_result```, and is now concrete.


### Implementation

 * Streamlined compute table implemenation using templates.

