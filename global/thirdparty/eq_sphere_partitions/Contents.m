% Recursive Zonal Equal Area Sphere Partitioning Toolbox.
% Release 1.10 2005-06-26
%
%Functions by category
%=====================
%
%Installation
%------------
%
%  install_eq_toolbox     Install toolbox
%  uninstall_eq_toolbox   Uninstall toolbox
%
%Recursive zonal equal area sphere partitions
%--------------------------------------------
%
% Partitions
%  eq_caps                Partition a sphere into to nested spherical caps
%  eq_regions             Recursive zonal equal area (EQ) partition of sphere
%
% Point sets
%  eq_point_set           Center points of regions of EQ partition, 
%                             in Cartesian coordinates
%  eq_point_set_polar     Center points of regions of an EQ partition, 
%                             in spherical polar coordinates
%
% Partition options
%  partition_options      Options for EQ partition
%
%Properties of recursive zonal equal area sphere partitions
%----------------------------------------------------------
%
% Diameter
%  eq_diam_bound          Maximum per-region diameter bound of EQ partition
%  eq_vertex_diam         Maximum vertex diameter of EQ partition
%
%  eq_diam_coeff          Coefficients of diameter bound and vertex diameter of
%                             EQ partition
%  eq_vertex_diam_coeff   Coefficient of maximum vertex diameter of EQ partition
%
%Hook for user-defined properties
%  eq_regions_property    Property of regions of an EQ partition
%
%Properties of EQ point sets
%---------------------------
%
% Energy and minimum distance
%  eq_energy_dist         Energy and minimum distance of an EQ point set
%  point_set_energy_dist  Energy and minimum distance of a point set
%
% Energy
%  eq_energy_coeff        Coefficient in expansion of energy of an EQ point set
%  point_set_energy_coeff Coefficient in expansion of energy of a point set
%  calc_energy_coeff      Coefficient of second term in expansion of energy
%
% Minimum distance
%  eq_min_dist            Minimum distance between center points of an 
%                             EQ partition
%  point_set_min_dist     Minimum distance between points of a point set
%
%  eq_dist_coeff          Coefficient of minimum distance of an EQ point set
%  point_set_dist_coeff   Coefficient of minimum distance of a point set
%  calc_dist_coeff        Coefficient of minimum distance
%
% Spherical cap packing density
%  eq_packing_density     Density of packing given by minimum distance of 
%                             EQ point set
%  point_set_packing_density  Density of packing given by minimum distance of 
%                             a point set
%  calc_packing_density   Density of packing given by minimum distance
%
% Hook for user-defined properties
%  eq_point_set_property  Property of an EQ point set
%
%Illustrations
%-------------
%
% Illustration of EQ Partitions of S^2 or S^3
%  show_s2_partition      3D illustration of an EQ partition of S^2
%  project_s2_partition   Use projection to illustrate an EQ partition of S^2
%  project_s3_partition   Use projection to illustrate an EQ partition of S^3
%
% Illustration of point sets on S^2 or S^3
%  show_r3_point_set      3D illustration of a point set
%  project_point_set      Use projection to illustrate a point set of 
%                             S^2 or S^3
% Illustration of a recursive zonal equal area sphere partition algorithm
%  illustrate_eq_algorithm  Illustrate the EQ partition algorithm
%
% Illustration options
%  illustration_options   Options for illustrations of EQ partitions
%
%Tests
%-----
%
%  eq_area_error          Total area error and max area error per region of an
%                             EQ partition
%
%Utilities
%---------
%
%  area_of_cap            Area of spherical cap
%  area_of_collar         Area of spherical collar
%  area_of_ideal_region   Area of one region of an EQ partition
%  area_of_sphere         Area of sphere
%  cart2polar2            Convert Cartesian to spherical polar coordinates on S^2
%  euc2sph_dist           Convert Euclidean to spherical distance
%  euclidean_dist         Euclidean distance between two points
%  fatcurve               Create a parameterized cylindrical surface
%  haslight               Check if axis handle has a light attached
%  ideal_collar_angle     Ideal angle for spherical collars of an EQ partition
%  illustration_options   Options for illustrations of EQ partitions
%  partition_options      Options for EQ partition
%  polar2cart             Convert spherical polar to Cartesian coordinates
%  sph2euc_dist           Convert spherical to Euclidean distance
%  spherical_dist         Spherical distance between two points on the sphere
%  sradius_of_cap         Spherical radius of spherical cap of given area
%  volume_of_ball         Volume of the unit ball
%
% If the previous text scrolled off the screen, try
% more on; help eq_sphere_partitions; more off;

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-26 $
% Function changed name from e2s to euc2sph_dist
% Function changed name from s2e to sph2euc_dist
% Function changed name from s2x to polar2cart
% Function changed name from x2s2 to cart2polar2
% Add new function fatcurve
% Add new function haslight
% Clean up descriptions
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.
