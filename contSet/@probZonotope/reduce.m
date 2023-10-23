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
%    Z1 = [1 2  4 2  5 3  4 2 1 3; ...
%         -3 2 -1 4 -3 2 -3 4 2 1];
%    Z2 = [3 -2 1 4 2; -1 1 0 2 1];
%    probZ = probZonotope(Z1,Z2);
%    probZred = reduce(probZ,'girard',3);
%
%    figure; hold on;
%    plot(probZ); plot(probZred);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       27-September-2007 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{probZ,'att','probZonotope'};
                {option,'str',{'girard','althoff'}};
                {order,'att','numeric','nonnan'}});

%reduce uncertain mean
Zred = reduce(zonotope(probZ.Z),option,order);
probZ.Z = Zred.Z;

if probZ.gauss~=1
    %reduce probabilistic part
    probZred=probReduce(probZ);
else
    probZred=probZ;
end

% ------------------------------ END OF CODE ------------------------------
