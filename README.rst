#########
SeeK-path
#########

``SeeK-path`` is a python module to obtain and visualize band paths in the
Brillouin zone of crystal structures. 

The definition of k-point labels follows crystallographic convention, as defined
and discussed in the `HPKOT paper`_. Moreover, the Bravais lattice is detected
properly using the spacegroup symmetry. Also the suggested band path provided
in the `HPKOT paper`_ is returned.
Systems without time-reversal and inversion-symmetry are also properly 
taken into account.

.. contents::

.. section-numbering::

===========
How to cite
===========
If you use this tool, please cite the following work:

- Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, *Band structure diagram 
  paths based on crystallography*, Comp. Mat. Sci. 128, 140 (2017)
  (`JOURNAL LINK`_, `arXiv link`_).
- You should also cite `spglib`_ that is an essential library used in the 
  implementation.

==============
How to install
==============
To install, use ``pip install seekpath``. It works both in python 2.7 and 
in python 3.5.

==========
How to use
==========
The main interface of the code is the python function 

    seekpath.get_path(structure, with_time_reversal, recipe, threshold)

You need to pass a crystal structure, a boolean flag (``with_time_reversal``) to say if time-reversal symmetry is present or not, and optionally, a recipe (currently only the string "HPKOT" is supported) and a numerical threshold.

The format of the structure is described in the function docstring. In particular,
It should be a tuple in the format 

``(cell, positions, numbers)``

where (if ``N`` is the number of atoms): 

- ``cell`` is a ``3x3`` list of floats (``cell[0]`` is the first lattice vector, ...); 
- ``positions`` is a ``Nx3`` list of floats with the atomic coordinates in scaled coordinates (i.e., w.r.t. the cell vectors);
- ``numbers`` is a length-``N`` list with integers identifying uniquely the atoms in the cell.

The output of the function is a dictionary containing, among other quantities, the k-vector coefficients, the suggested band path, whether the system has inversion symmetry, the crystallographic primitive lattice, the reciprocal primitive lattice.
A detailed description of all output information and their format can be found in the function docstring.

---------------------------------------------------------------
A warning on how to use (and crystal structure standardization)
---------------------------------------------------------------
SeeK-path standardizes the crystal structure 
(e.g., rotates the tetragonal system so that the *c* axis is along *z*, 
etc.) and can compute the suggested band paths only of standardized 
(crystallographic) primitive cells. Therefore, the 
**correct approach to use this tool is the following**:

1. You first find the standardized primitive cell with SeeK-path (returned in
   output) and store it somewhere, together with the k-point coordinates
   and suggested band path

2. You then run all your calculations using the standardized primitive cell

If you already have done calculations with a non-standardized cell, you will
then need to figure out how to remap the labeled k-points in the choice of
cell you did.

Note: SeeK-path returns, in output, a transformation matrix. You can use
this to understand the relationship between your input cell and the one
standardized by SeeK-path. More details can be found below in the section
*Transformation matrix*.

---------------
Explicit k path
---------------

You might also be interested in the function 

     seekpath.get_explicit_k_path

that has a very similar interface, that produces an explicit list of k-points along
the suggested band path. The function has the same interface as ``get_path``, but 
has also an additional optional parameter ``reference_distance``, that is used as a reference target distance between neighboring k-points along the path. More detailed information can be found in the docstrings.

---------------------
Transformation matrix
---------------------
Before providing the suggested band path, SeeK-path will standardise the crystal structure (using routines in spglib, except for the case of triclinic cells).

The standardisation consists of two steps:

1. **rearrangement of the lattice vectors as a linear combination of them** (this, e.g., involves permutations, inversion of axes, changing a vector with a linear combination of them, ...). For instance, this will reshuffle the three orthogonal vectors of a tetragonal cell, so that the third vector is the one with different length. Or, for a supercell, it will rescale the cell to get the conventional cell.

2. **rotation of Cartesian axes**: vectors, their lengths and relative angles are now chosen. The triplet of vectors can now be rotated in Cartesian space, e.g., to have the first vector along the ``x`` axis, the second along ``y`` etc. Note that, in this step, if the cell is refined, cell lengths and relative angles can be slightly adjusted (e.g. the length of the three vectors is set to be the same for almost-cubic systems, and the angles to be exactly 90 degrees even if they weren't so).

The ``get_path`` function of SeeK-path returns, on output, a ``conv_transformation_matrix`` ``T``, that contains information on step 1. 

This is defined as follows::

  cell_orig = T * cell_std

where ``cell_orig`` is the original cell provided in input by the user, ``cell_std`` is the *standardized, conventional cell* and ``*`` represents matrix multiplication. Note that the cells are written in the same format as in input to SeeK-path: the first row ``cell[0,:]`` is the first cell vector, etc.

* *Note 1*: that this transformation matrix relates the original cell with the *conventional* cell of the crystal, rather than the primitive. SeeK-path provides also the transformation from conventional to primitive independently, that allows you to get the total transformation matrix.

* *Note 2*: that this transformation matrix is the transpose of the one returned by spglib.

The matrix ``T`` can be used to know how the original cell is related to the final one (only for what concerns step 1).

For step 2, you can calculate the rotation matrix ``R`` multiplying the inverse of the result of the formula above containing the transformation matrix (using the output cell of SeeK-path), by the actual input to SeeK-path. The result ``R`` is normally a rotation matrix. However, special care must be used if the structure has been refined, because in this case the matrix obtained with this approach will not be exactly orthogonal.

This means that, globally, the ``cell_input`` given as input to SeeK-path 
is given by::

  cell_input = T * cell_std * R

where ``cell_std`` is the *conventional cell* output by SeeK-path, and ``R``
is a rotation only if ``cell_input`` was already refined.




=================
AiiDA integration
=================

If you use AiiDA (www.aiida.net), you might be interested in replacing the above
functions with the following wrappers, instead:

    seekpath.aiidawrappers.get_path 
    
    seekpath.aiidawrappers.get_explicit_k_path 

The function interfaces are very similar, but the advantage is that these functions expect an AiiDA structure as input (instead of a tuple) and return AiiDA structures and KpointsData classes instead of lists and tuples, where appropriate.
Also in this case, additional information is found in the docstrings.


=======
License
=======

The code is open-source (licensed with a MIT license, see LICENSE.txt).

===================
Online service/tool
===================

In this repository we also provide the code to deploy a online service for 
the visualization of the band paths and primitive cells of the crystal 
structures. A live demo is currently hosted on the `MaterialsCloud`_ web portal.

The following is a screenshot of the selection window:

.. image:: https://raw.githubusercontent.com/giovannipizzi/seekpath/master/webservice/screenshots/selector.png
     :alt: SeeK-path web service selection window
     :width: 50%
     :align: center

And the following is a screenshot of the main output window, showing the Brillouin zone, the primitive crystal structure, the coordinates of the k-points and the suggested band path.

.. image:: https://raw.githubusercontent.com/giovannipizzi/seekpath/master/webservice/screenshots/mainwindow.png
     :alt: SeeK-path web service main output
     :width: 50%
     :align: center

.. _HPKOT paper: http://dx.doi.org/10.1016/j.commatsci.2016.10.015
.. _JOURNAL LINK: http://dx.doi.org/10.1016/j.commatsci.2016.10.015
.. _arXiv link: https://arxiv.org/abs/1602.06402
.. _spglib: http://atztogo.github.io/spglib/
.. _MaterialsCloud: http://www.materialscloud.org/tools/seekpath/
