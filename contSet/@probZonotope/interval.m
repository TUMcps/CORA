function Int = interval(probZ)
% interval - Overapproximates a probabilistic zonotope by an interval
%
% Syntax:
%    Int = interval(probZ)
%
% Inputs:
%    probZ - probabilistic zonotope object
%
% Outputs:
%    Int - interval object
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2,2);
%    Int = interval(probZ);
%
% Other m-files required: interval (constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       14-September-2006 
% Last update:   22-March-2007
%                25-July-2016 (intervalhull replaced by interval)
%                27-August-2019 (generators, center)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%extract center and generators
c = center(probZ);
G = generators(probZ);

%determine left and right limit
leftLimit = c - sum(abs(G),2);
rightLimit = c + sum(abs(G),2);

%instantiate interval
Int = interval(leftLimit,rightLimit);

% ------------------------------ END OF CODE ------------------------------
