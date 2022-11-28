function p = randPoint(probZ,varargin)
% randPoint - generates a random point of a probabilistic 
% zonotope
%
% Syntax:  
%    p = randPoint(probZ)
%    p = randPoint(probZ,N)
%    p = randPoint(probZ,N,type)
%    p = randPoint(probZ,'all','extreme')
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
% See also: zonotope/randPoint

% Author:       Matthias Althoff
% Written:      17-July-2020
% Last update:  19-August-2022 (MW, integrate standardized pre-processing)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_randPoint('probZonotope',probZ,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of vars
    p = vars{1}; return
else
    % assign variables
    probZ = vars{1}; N = vars{2}; type = vars{3}; pr = vars{4};
    % sampling with 'gaussian' is done in contSet method
    if strcmp(type,'gaussian')
        p = randPoint@contSet(probZ,N,type,pr); return
    end
end

% generate random points within the corresponding zonotope 
c = randPoint(zonotope(probZ.Z),N,type);

% generate random point from normal distribution
p = zeros(size(c));
for i = 1:size(p,2)
    p(:,i) = mvnrnd(c(:,i),sigma(probZ))';
end

%------------- END OF CODE --------------