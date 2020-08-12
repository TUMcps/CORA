function probZred = reduce(probZ,option,order)
% reduce - Reduces the order of a probabilistic zonotope
%
% Syntax:  
%    probZred = reduce(probZ,option,order)
%
% Inputs:
%    probZ - probabilistic zonotope object
%    option - available options:
%             - 'girard'
%             - 'althoff' (order set to 1)
%    order - order of reduced zonotope
%
% Outputs:
%    probZred - reduced zonotope
%
% Example:
%    Z1 = rand(2,20);
%    Z2 = -1+2*rand(2,5);
%    probZ = probZonotope(Z1,Z2);
%    probZred = reduce(probZ,'girard',3);
%
%    plot(probZ,'dark'); hold on;
%    plot(probZred,'light');
%
% Other m-files required: none
% Subfunctions: reduceGirard, reduceAlthoff
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       27-September-2007 
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

%reduce uncertain mean
Zred = reduce(zonotope(probZ.Z),option,order);
probZ.Z = Zred.Z;

if probZ.gauss~=1
    %reduce probabilistic part
    probZred=probReduce(probZ);
else
    probZred=probZ;
end

%------------- END OF CODE --------------