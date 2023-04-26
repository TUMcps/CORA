Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox.
Release 1.10 2005-06-26
Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
README.txt

For licensing, see COPYING
For references, see AUTHORS
For revision history, see CHANGELOG.

NOTE:
This file has lines terminated by CR-LF for use with DOS and Windows Notepad.
On Windows, to read AUTHORS, COPYING and CHANGELOG, use another editor, such as 
Wordpad or the Matlab editor.

>> What is the Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox?

The Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox is a suite of
Matlab functions. These functions are intended for use in exploring different
aspects of EQ sphere partitioning.

The functions are grouped into the following groups of tasks:
1. Create EQ partitions
2. Find properties of EQ partitions
3. Find properties of EQ point sets
4. Produce illustrations
5. Test the toolbox
6. Perform some utility function

>> What is an EQ partition?

An EQ partition is a partition of S^dim [the unit sphere in the dim+1 Euclidean
space R^(dim+1)] into a finite number of regions of equal area. The area of
each region is defined using the Lebesgue measure inherited from R^(dim+1).

The diameter of a region is the sup of the Euclidean distance between any two
points of the region. The regions of an EQ partition have been proven to have 
small diameter, in the sense that there exists a constant C(dim) such that the
maximum diameter of the regions of an N region EQ partition of S^dim is bounded
above by C(dim)*N^(-1/dim).

>> What is an EQ point set?

An EQ point set is the set of center points of the regions of an EQ partition.
Each region is defined as a product of intervals in spherical polar coordinates.
The center point of a region is defined via the center points of each interval,
with the exception of spherical caps and their descendants, where the center
point is defined using the center of the spherical cap.

>> What's new?

Release 1.00 2005-02-13 was the first version of the EQ sphere partitioning
code to be packaged in toolbox form.
Release 1.10 2005-06-26 is the first public release of the code.
For details of changes from earlier versions of the EQ partitioning code, 
and changes from 1.00 to the current release, see CHANGELOG.

>> Which versions of Matlab can I use?

This toolbox has been tested with Matlab versions 6.5 and 7.0.1 on Linux,
and 6.5.1 on Windows.

>> How do I install the Recursive Zonal Equal Area Sphere Partitioning Toolbox?

This toolbox is organized into a number of directories. To use it effectively,
these directories need to be on your Matlab path every time you start Matlab.
You will therefore need to install the toolbox before using it.

To do this,
1. Unzip the file eqsp-1.10.zip into the directory where you
   want the toolbox to reside. This will create the subdirectory
   eq_sphere_partitions.
2. Run Matlab, change directory to eq_sphere_partitions and then run
   install_eq_toolbox.
For more information, see INSTALL.txt.

>> What documentation is available?

The user documentation includes the help comments in each Matlab M file,
plus the extra files Contents.m, AUTHORS, CHANGELOG, COPYING, INSTALL.txt
and README.txt.

The extra files contain the following:
Contents.m:    A brief description of those functions in the toolbox which are
               visible to end users.
AUTHORS:       Authors, acknowledgements and references.
CHANGELOG:     Revision history.
COPYING:       Software license terms.
INSTALL.txt:   Installation instructions in DOS CR-LF text format.
README.txt:    This file.

>> How do I get help?

To see a brief description of the functions in the toolbox, enter the command
HELP EQ_SPHERE_PARTITIONS (in lower case).
The command WHAT lists all the functions in your current directory.

For each function, the command HELP FUNCTION, where FUNCTION is
the name of the function, will give the help for the function.

The help format is:
1. Summary line.
2. Syntax of the function call.
3. Description, usually with a description of each argument.
4. Notes.
5. Example of use of the function.
6. See also.

For example, here is the result of HELP EQ_POINT_SET:
-------------------------------------------------------------------------------
 EQ_POINT_SET Center points of regions of EQ partition, in Cartesian coordinates

 Syntax
  points_x = eq_point_set(dim,N,options);

 Description
  POINTS_X = EQ_POINT_SET(dim,N) does the following:
  1) uses the recursive zonal equal area sphere partitioning algorithm to
  partition S^dim (the unit sphere in dim+1 dimensional space) into N regions
  of equal area and small diameter, and
  2) sets POINTS_X to be an array of size (dim+1 by N), containing the center
  points of each region.
  Each column of POINTS_X represents a point of S^dim, in Cartesian coordinates.

  The arguments dim and N must be positive integers.

  POINTS_X = EQ_POINT_SET(dim,N,'offset','extra') uses experimental extra offsets
  for S^2 and S^3 to try to minimize energy.

  POINTS_X = EQ_POINT_SET(dim,N,extra_offset) uses experimental extra offsets if
  extra_offset is true or non-zero.

 Notes
  Each region is defined as a product of intervals in spherical polar
  coordinates. The center point of a region is defined via the center points
  of each interval, with the exception of spherical caps and their descendants,
  where the center point is defined using the center of the spherical cap.

  If dim > 3, extra offsets are not used.
  For more details on options, see help partition_options.

 Examples
  > points_x = eq_point_set(2,4)
  points_x =
           0    0.0000   -0.0000    0.0000
           0    1.0000   -1.0000         0
      1.0000    0.0000    0.0000   -1.0000

  > size(points_x)
  ans =
       3     4

 See also
  PARTITION_OPTIONS, EQ_POINT_SET_POLAR, EQ_REGIONS, S2X
-------------------------------------------------------------------------------

>> What is the input? What is the output?

The help for each function briefly describes the input and the output, as per
the example for partsphere, above.

>> Which file to begin with?

You need to find a function which does what you want to do. Examples:

1. Create EQ partitions

o  Create an array in Cartesian coordinates representing the `center' points
of an EQ partition of S^dim into N regions:

> points_x = eq_point_set(dim,N);

o  Create  an array in spherical polar coordinates representing the `center'
points of an EQ partition of S^dim into N regions:

> points_s = (eq_point_set_polar(dim,N);

o  Create an array in polar coordinates representing the regions of an EQ
partition of S^dim into N regions:

> regions = eq_regions(dim,N);

2. Find properties of EQ partitions

o Find the (per-partition) maximum diameter bound of the EQ partition of S^dim
into N regions:

> diam_bound = eq_diam_bound(dim,N);

3. Find properties of EQ point sets

o Find the r^(-s) energy and min distance of the EQ `center' point sets of
S^dim for N points:

> [energy,dist] = eq_energy_dist(dim,N,s);

4. Produce an illustration

o Use projection to illustrate the EQ partition of S^2 into N regions:

> project_s2_partition(N);

>> Is the toolbox for use with S^2 and S^3 only? What is the maximum dimension?

In principle, any function which has dim as a parameter will work for any
integer dim >= 1. In practice, for large d, the function may be slow, 
or may consume large amounts of memory.

>> What is the range of the number of points, N?

In principle, any function which takes N as an argument will work with any
positive integer value of N. In practice, for very large N, the function may
be slow, or may consume large amounts of memory.

>> What are the options for visualizing points or equal area regions?

A number of different illustrations are available:

o Use a 3D plot to illustrate the EQ partition of S^2 into N regions:
> show_s2_partition(N);

o Use projection to illustrate the EQ partition of S^2 into N regions:
> project_s2_partition(N);

o Use projection to illustrate the EQ partition of S^3 into N regions.
> project_s3_partition(N);

o Illustrate the EQ algorithm for the partition of S^dim into N regions.
> illustrate_eq_algorithm(dim,N);

See the help for these functions for details.
