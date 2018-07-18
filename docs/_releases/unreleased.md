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

The compute table interface has changed significantly.
This will affect anyone who implements their own operations.
Essentially, the compute table is now much more aware of the contents
of each compute table entry (in terms of the types of items in an entry).
The long term goals of this change include
simplification (and better flexibility flexibility) of operation code
and eventually to allow for efficient compression of compute table entries.
   

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

  * Instead of returning the result, the ```compute_table::find``` 
    method now expects the result to be passed as an argument
    and will be filled in.
    Operations now pre-allocate results for this purpose
    (only one result is needed per entry type).

  * Removed ```OperationMap``` as a compute table option.


### Implementation

 * Streamlined compute table implemenation using templates.

