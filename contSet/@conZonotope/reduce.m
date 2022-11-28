function cZ = reduce(cZ,method,order,varargin)
% reduce - reduce the number of generators of a constrained zonotope object
%
% Syntax:  
%    cZ = reduce(cZ,method,order)
%    cZ = reduce(cZ,method,order,redOptions)
%
% Inputs:
%    cZ - conZonotope object
%    method - zonotope reduction method (i.e 'girard, 'combastel', etc.). 
%    order - desired degree-of-freedom order
%    redOptions - additional settings for the zonotope reduction method 
%                 (i.e. filterLength, alg, etc.)
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    Z = [0 1 0 1 1;0 1 2 -1 0];
%    A = [-2 1 -1 2]; b = 2;
%    cZ = conZonotope(Z,A,b);
%    cZred = reduce(cZ,'girard',1);
%
%    figure; hold on;
%    plot(cZ,[1,2],'b')
%    plot(cZred,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce, reduceConstraints
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:      11-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% Reduce the number of generators. Implementation according to Sec. 4.3
% in reference paper [1]

% parse input arguments
redOptions = {[]};
if nargin >= 4
    redOptions = varargin(1:end);
end

% remove the trivial constraints [0 0 ... 0]*beta = 0
cZ = deleteZeros(cZ);

% object properties
A = cZ.A;
b = cZ.b;
m = size(cZ.Z,1);

if ~isempty(A)

    % up-lift a c-zonotope to an ordinary zonotope
    Z_up = zonotope([cZ.Z; -b, A]);

    % calculate the reduced order for the lift-up-zonotope that is 
    % required to obtain the desired degree-of-freedom order
    nc = size(cZ.A,1);
    order = max(1,(order*m + nc)/(m+nc));

    % reduce the lift-up zonotope with convential zonotope reduction 
    % technique
    Z_up = reduce(Z_up,method,order,redOptions{:});

    % down-lift to a constrained zonotope
    C = Z_up.Z;
    cZ.Z = C(1:m, :);
    cZ.A = C(m+1:end, 2:end);
    cZ.b = -C(m+1:end, 1);

else

    % reduce the lift-up zonotope with convential zonotope reduction 
    % technique
    zRed = reduce(zonotope(cZ.Z),method,order,redOptions{:});
    cZ.Z = zRed.Z;   
end

%------------- END OF CODE --------------