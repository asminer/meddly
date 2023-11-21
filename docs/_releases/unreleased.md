---
title: New since last release
shorttitle: Since last release
number: changes
layout: single
---

### Interface Changes

* Moved ```node_headers``` inner classes to ```arrays.h```

* Created ```node_marker``` object, for mark and sweep
    (or just various mark applications).

* Created ```dot_maker``` object in file ```io_dot.h```,
    for creating dot files from MDDs.

* Methods ```expert_forest::writeNodeGraphPicture```
    and ```dd_edge::writePicture```
    are now deprecated; use object ```dot_maker``` instead.

### Implementation

