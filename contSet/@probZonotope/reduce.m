function [pZred]=reduce(pZ,option,order)
% reduce - Reduces the order of a probabilistic zonotope
% option 'girard': Use order reduction technique by Antoine Girard
% option 'althoff': Use order reduction technique by Matthias Althoff where
% order is set to 1
%
% Syntax:  
%    [Zred]=reduce(Z,option,order)
%
% Inputs:
%    pZ - probabilistic zonotope object
%    option - 'girard' or 'althoff'
%    order - order of reduced zonotope
%
% Outputs:
%    pZred - reduced zonotope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: reduceGirard, reduceAlthoff
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 27-September-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%reduce uncertain mean
Zred=reduce(zonotope(pZ.Z),option,order);
pZ.Z=get(Zred,'Z');

if pZ.gauss~=1
    %reduce probabilistic part
    pZred=probReduce(pZ);
else
    pZred=pZ;
end

%------------- END OF CODE --------------