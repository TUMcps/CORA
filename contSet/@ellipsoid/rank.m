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

if representsa_(E,'emptySet',eps)
    r = 0;
else
    r = zeros(size(E));
    for i=1:numel(E)
        % Find minimum svd threshold using reciprocal
        % condition number
        d = svd(E(i).Q);
        mev_th = d(1)*E(i).TOL;
        r(i) = sum(d>0 & d>=mev_th);
    end
end

% ------------------------------ END OF CODE ------------------------------
