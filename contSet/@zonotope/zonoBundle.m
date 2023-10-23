function zB = zonoBundle(Z)
% zonoBundle - convert a zonotope object into a zonotope bundle object
%
% Syntax:
%    zB = zonoBundle(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    zB - zonoBundle object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope/zonoBundle

% Authors:       Niklas Kochdumper
% Written:       26-November-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
zB = zonoBundle({Z});

% ------------------------------ END OF CODE ------------------------------
