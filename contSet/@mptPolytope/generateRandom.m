function poly = generateRandom(varargin)
% generateRandom - Generates a random zonotope
%
% Syntax:  
%    poly = generateRandom()
%    poly = generateRandom(dim)
%    poly = generateRandom(dim,cen)
%
% Inputs:
%    dim      - dimension
%    cen      - (optional) center
%
% Outputs:
%    poly - random mptPolytope
%
% Example: 
%    poly = mptPolytope.generateRandom(2);
%
%    plot(poly);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/generateRandom

% Author:       Niklas Kochdumper
% Written:      05-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments 
    n = rand(1,5);

    if nargin >= 1 && ~isempty(varargin{1})
       n = varargin{1}; 
    end

    cen = -10 + 20*rand(n,1);

    if nargin >= 2 && ~isempty(varargin{2})
       cen = varargin{2}; 
    end

    % generate random normal distribution
    A = rand(n,n);

    sigma = 0.5*(A+A') + n*eye(n);
    points = mvnrnd(cen,sigma,n + 100)';

    % generate random mptPolytope
    poly = mptPolytope.enclosePoints(points);

%------------- END OF CODE --------------