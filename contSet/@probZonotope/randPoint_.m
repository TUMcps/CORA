function p = randPoint_(probZ,N,type,varargin)
% randPoint_ - generates a random point of a probabilistic zonotope
%
% Syntax:
%    p = randPoint_(probZ)
%    p = randPoint_(probZ,N)
%    p = randPoint_(probZ,N,type)
%    p = randPoint_(probZ,'all','extreme')
%
% Inputs:
%    probZ - probabilistic zonotope object
%    N - number of random points
%    type - type of the random point ('standard' or 'extreme')
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    Z1 = [10 1 -2; 0 1 1];
%    Z2 = [0.6 1.2; 0.6 -1.2];
%    probZ = probZonotope(Z1,Z2);
%    p = randPoint(probZ);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint, zonotope/randPoint_

% Authors:       Matthias Althoff
% Written:       17-July-2020
% Last update:   19-August-2022 (MW, integrate standardized pre-processing)
% Last revision: 27-March-2023 (MW, rename randPoint_)

% ------------------------------ BEGIN CODE -------------------------------

% generate random points within the corresponding zonotope 
c = randPoint_(zonotope(probZ.Z),N,type);

% generate random point from normal distribution
p = zeros(size(c));
for i = 1:size(p,2)
    p(:,i) = mvnrnd(c(:,i),sigma(probZ))';
end

% ------------------------------ END OF CODE ------------------------------
