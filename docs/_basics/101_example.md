---
title: Build MDD from function
shorttitle: Build MDD from function
mathjax: true
---

Build an MDD for Boolean function $f = x \vee (y \wedge z)$

## <span style="color:blue">  Truth Table </span>

| $x$ | $y$ | $z$ | $f$ |	
|--|--|--|--|
|0|0|0|0|
|0|0|1|0|
|0|1|0|0|
|0|1|1|1|
|1|0|0|1|
|1|0|1|1|
|1|1|0|1|
|1|1|1|1|


## <span style="color:blue"> Preliminaries <span>

- An MDD in MEDDLY, is characterized by an object of class forest. Each forest is built upon a domain of variables and belongs to certain category which is specified as per the set/function that needs to be represented.

 - To start with building MDD, each $v$ of the variables from the function is assigned a level, $k$ such that $var(k) = v$. 
 In this example, $x = var(3), y = var(2), z = var(1).$  

- Since the function here is Boolean, the domain of each variable is [0,1]. The output of the function is also Boolean, which is defined by the terminals of the MDD.

### <span style="color:blue"> Create Domain <span>

```c
domain* d = createDomain();
d->createVariablesBottomUp(var_bounds, num_of_vars);
```

> num_of_vars = 3
> var_bounds : var_bounds[$i$] = 2, $\forall i\in$ [0 , num_of_vars-1]

### <span style="color:blue"> Create Forest on the Domain <span>

```
forest* mdd = d->createForest(false, forest::BOOLEAN, forest::MULTI_TERMINAL);
```

> createForest() takes parameter for whether mdd is for relations, range of terminals and edge labelling range.

### <span style="color:blue"> Create an MDD edge from an element of truth table </span>

```c
dd_edge first(mdd);
mdd->createEdge(element, 1, first);
```

> element : element[$i$] = $var(i)$.val, $\forall i \in$ [1,3]. 
> Default value of element[0] (terminal) is true.

###  <span style="color:blue"> Add multiple elements into MDD. </span>
- Via creatEdge() only : User can pass all elements and MEDDLY will take care of  adding them to the forest.
 All elements can be stored in an array and passed into createEdge(), like above example :
 
```c
dd_edge all(mdd);
mdd->createEdge(elementList, 8, all);
```
>  elementList : elementList[$j$] = element,  $\forall j \in$ [1,8]

- Via createEdge() and UNION operation : User can pass elements one-by-one and explicitly union them.
 
 ```c
dd_edge result(mdd);
for each element from truth-table : 
    mdd->createEdge(element,1,first);
    apply(UNION, first, result, result);
```

> result.show(std::cout, 0) will display the stored mdd.
