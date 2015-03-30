Selector command type
======================

The _selector command type_ defines selector objects in CPPPO. At the moment there are three _command specifiers_ for this _command type_:

* _cellIJK_ 

* _cellUnstruct_ 

* _filter_

A command with the _cellIJK_ specifier represents a cell selector based on IJK indexing. A command with the _cellUnstruct_ specifier represents a cell selector based on unstructured mesh. There are significant differences and
usages between the twos:

* A _cellIJK_ selector is based on the assumption that all the cells are equal (same volumes, same lengths) and hexahedrals. Furthermore, the domain is hexahedral itself. If this is the case, any cell can be identifyed using an IJK notation. This results in a fast algorithms because, if a particular region in the domain is given, the number and ids of the cells inside that region are calculated using simple relations.

* A _cellUnstruct_ selector does not pose any limitation on the shape and lenghts of cells and domain. The only drawback of this selector is the speed which is lower compared to the _cellIJK_ selector when the filtersize is small. 

Both these selectors require the domain to be properly decomposed. Sub-domains can have any shape (if the _cellUnstruct_ seector is used) but the hexahedrals built using maxima and minima of every subdomain should not overlap. In other words, processor-processor boundaries should be perpendicular or parallel to the main axes. 

A command with the _filter specifier_ represents a _filter_. _filters_ have to be defined in the c3po.json file. For filters, cells are selected once their center is within the filter distance (specified in each direction, see below) from the master cell, i.e., the cell for which the filtered value is calculated. This means, that the filter width specified below is half of the box size over which the filtering (i.e., spatial averaging) is performed.

WARNING: It is recommended to define the filter width such that it extends a certain distance further than the cell center. In such a way, problems due to not selecting the cell can be avoided. Anyhow, CPPPO will add a small value (i.e., the value for `filterWidthTolerance` specified in the `mesh.json` file) to the filter width in order ensure the cell is selected in case the filter matches exactly the distance between cell centers.

Syntax
------
```
...

 "myFilter0":	  { 
                     "CoordSys": 0, 
                     "x":  1.0 , 
                     "y":  1.0 ,
                     "z":  1.0 
                  },
...

```

Every _filter_ requires the following entry:

* "CoordSys": requires an integer value. If set to 0 the cell selection will be performed in cartesian coordinates and if set to 1, spherical coordinates will be used instead.

In the case of a cartesian coordinate system, three more entries are required:

*  "x": requires a double value. This entry specifies the spatial extension of the filter into the x direction.

*  "y": requires a double value. This entry specifies the spatial extension of the filter into the y direction.

*  "z": requires a double value. This entry specifies the spatial extension of the filter into the z direction.

In the case of a spherical coordinate system, one more entry is required:

*  "r": requires a double value. This entry specifies the spatial extension of the filter into the radial direction.

For every _filter_ it is necessary to have a corresponding _cellIJK_ (or _cellUnstruct_) specifier. The corresponding _cellIJK_ (or _cellUnstruct_) specifier is selected based on the ordering in the `c3po.input` script. 



Example
-------

In the `c3po.input` file:

```
...

selector     cellUnstruct    cellSelector0
selector     cellIJK         cellSelector1

selector     filter          myFilter0
selector     filter          myFilter1 

...

```
While in the `c3po.json` file:

```
...

 "myFilter0":	  { 
                     "CoordSys": 0, 
                     "x":  1.0 , 
                     "y":  1.0 ,
                     "z":  1.0 
                  },
                  
 "myFilter1":	  { 
                     "CoordSys": 0, 
                     "x":  2.0 , 
                     "y":  1.0 ,
                     "z":  3.0 
                  },
...
```
In the previous example, _myFilter0_ will use _cellSelector0_ (the unstructured selector) and _myFilter1_ will use _cellSelector1_ (the IJK selector).
Note that the choice of the selector only depends from the mesh and the `c3po.input` should be:

```
...

selector     cellIJK         cellSelector0
selector     cellIJK         cellSelector1

selector     filter          myFilter0
selector     filter          myFilter1 

...

```         
In the case of a structured mesh without any refinement or a simple geometry. Alternatively it should have been:
```
...

selector     cellUnstruct         cellSelector0
selector     cellUnstruct         cellSelector1

selector     filter          myFilter0
selector     filter          myFilter1 

...

```      
In the case of an unstructured mesh or a complex geometry.

Go back
-----------
 - [commands](10_commandTypes.md) 
