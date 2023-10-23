function Z = box(Z)
% box - computes an enclosing axis-aligned box; the result is equivalent to
%    a conversion to intervals but yields a zonotope representation
%
% Syntax:
%    Z = box(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z = zonotope([1;-1],[-3 2 1; -1 0 3]);
%    B = box(Z);
%    
%    figure; hold on;
%    plot(Z);
%    plot(B,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       09-March-2009
% Last update:   27-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine new generator matrix
G = diag(sum(abs(Z.G),2));

% instantiate axis-aligned zonotope
Z.G = G;

% ------------------------------ END OF CODE ------------------------------
