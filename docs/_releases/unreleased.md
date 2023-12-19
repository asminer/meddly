---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Updates to unpacked nodes

The ```unpacked_node``` interface has been updated.
For efficiency, an unpacked node is now permanently attached to
one forest, and there is a free list of unpacked nodes for each forest.
Also, the interface has been updated to increase readability.


### Updates to forests

* All methods have been moved from ```expert_forest``` class
    to ```forest``` class.

* All instances of ```expert_forest``` have been replaced by ```forest```.


### Deprecated methods

Uncomment ```ALLOW_DEPRECATED_0_17_4``` in ```defines.h```
to use these deprecated methods.
This is intended to help developers migrate to the new interface.
Deprecated methods appear at the bottom of ```forest.h``` and ```unpacked.h```.

