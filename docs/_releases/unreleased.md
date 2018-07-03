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

### Expert Interface Changes

 * Class ```compute_table::search_key```
   has been renamed ```compute_table::entry_key```, and is now concrete.

 * Class ```compute_table::search_result```
   has been renamed ```compute_table::entry_result```, and is now concrete.

 * Class ```compute_table::entry_builder```
   has been removed.  A compute table entry is now added by giving
   the key portion and the result portion.

  * New class, ```compute_table::entry_type```
    with type information for compute table entries.

  * The list of recycled ```entry_keys``` is now maintained in 
    class ```compute_table```.
    When implementing an operation, old code:
    ```c++
    compute_table::entry_key* CTsrch = useCTkey();
    CTsrch->reset();
    // ...
    doneCTkey(CTsrch);
    ```
    should be replaced with:
    ```c++
    compute_table::entry_key* CTsrch = CT->useEntryKey(this);
    // ...
    CT->recycle(CTsrch);
    ```
  * Removed ```OperationMap``` as a compute table option.


### Implementation

 * Streamlined compute table implemenation using templates.

