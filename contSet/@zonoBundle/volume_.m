function vol = volume_(zB,varargin)
% volume_ - Computes the volume of a zonotope bundle
%
% Syntax:
%    vol = volume_(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    vol - volume
%
% Example: 
%    Z1 = zonotope([1;1], [1 1; -1 1]);
%    Z2 = zonotope([-1;1], [1 0; 0 1]);
%    zB = zonoBundle({Z1,Z2}); 
%    volume(zB)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/volume

% Authors:       Matthias Althoff
% Written:       02-February-2011 
% Last update:   18-August-2022 (MW, include standardized preprocessing)
% Last revision: 27-March-2023 (MW, rename volume_)

% ------------------------------ BEGIN CODE -------------------------------

%obtain polytope of zonotope bundle
P = polytope(zB);

%compute volume
vol = volume(P);

% ------------------------------ END OF CODE ------------------------------
