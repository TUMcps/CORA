function r = rank(E)
% rank - computes the dimension of the affine hull of an ellipsoid
%
% Syntax:
%    r = rank(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    r - dimension of the affine hull
%
% Example: 
%    E = ellipsoid([1,0;0,1]);
%    r = rank(E) 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       16-March-2021
% Last update:   04-July-2022 (VG, allow class array input)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
if representsa_(E,'emptySet',eps)
    r = 0;
    return
end

% Find minimum svd threshold using reciprocal condition number
d = svd(E.Q);
mev_th = d(1)*E.TOL;
r = sum(d > 0 & d >= mev_th);

% ------------------------------ END OF CODE ------------------------------
