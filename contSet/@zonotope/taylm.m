function t = taylm(Z)
% taylm - enclose a zonotope object with a Taylor model
%
% Syntax:
%    t = taylm(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    t - taylm object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Niklas Kochdumper
% Written:       13-August-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
% create taylor models for factors
m = size(Z.G,2);
dom = interval(-ones(m,1),ones(m,1));
t = taylm(dom);

% create taylor model for the zonotope
t = Z.c + Z.G*t;

% ------------------------------ END OF CODE ------------------------------
