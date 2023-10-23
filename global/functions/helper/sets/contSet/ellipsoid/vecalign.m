function T = vecalign(x,y)
% vecalign - computes T such that x || T*y (x parallel to T*y)
%
% Syntax:
%    T = vecalign(x,y)
%
% Inputs:
%    x - vector
%    y - vector
%
% Outputs:
%    T - matrix
%
% Example:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

x = x/norm(x);
y = y/norm(y);
[U1,~,~] = svd(x);
[U2,~,~] = svd(y);
T = U1*U2';

% ------------------------------ END OF CODE ------------------------------
