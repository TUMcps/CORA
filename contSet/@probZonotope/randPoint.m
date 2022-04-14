function p = randPoint(obj,varargin)
% randPoint - generates a random point of a probabilistic 
% zonotope
%
% Syntax:  
%    p = randPoint(obj)
%    p = randPoint(obj,N)
%    p = randPoint(obj,N,type)
%    p = randPoint(obj,'all','extreme')
%
% Inputs:
%    obj - probabilistic zonotope object
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
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    N = 1;
    type = 'standard';
    if nargin > 1 && ~isempty(varargin{1})
       N = varargin{1}; 
    end
    if nargin > 2 && ~isempty(varargin{2})
       type = varargin{2}; 
    end
    
    % generate random points within the corresponding zonotope 
    c = randPoint(zonotope(obj.Z),N,type);

    % generate random point from normal distribution
    p = zeros(size(c));
    for i = 1:size(p,2)
        p(:,i) = mvnrnd(c(:,i),sigma(obj))';
    end
end

%------------- END OF CODE --------------