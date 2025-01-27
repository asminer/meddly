---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Interface Changes

Several changes to the ```unpacked_node``` class:

* An ```unpacked_node``` object is now permanently tied to
    a forest, and to a storage method (full or sparse).

* The constructor is now private; new instances should be created
    using method ```New()```

* Methods ```newFull()``` and ```newSparse()``` are deprecated,
    in favor of ```newWritable()```

* For consistency, creating an unpacked node from a reduced
    node in a forest should be done via method ```initFromNode()```
    or ```newFromNode()```, instead of methods in class ```forest```.

### Implementation

