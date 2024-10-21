function c = center(SpS)
% center - Computes a point in a spectrahedron. Currently, this produces
%    only some feasible point, though this can be more accurate if SpS has
%    been constructed e.g., as a polytope. Better methods will be
%    developped in the future.
%
% Syntax:
%    c = center(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    c - Point in SpS
%
% Example:
%    A0 = eye(4);
%    A1 = blkdiag([1 0;0 -1],zeros(2));
%    A2 = blkdiag(zeros(2),[1 0;0 -1]);
%    SpS = spectraShadow([A0 A1 A2]);
%    c = center(SpS);
% 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       12-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if representsa_(SpS,'emptySet',1e-8)
    c = [];
    return
end

if ~isempty(SpS.center.val)
    c = SpS.center.val;
    return
end

c_spectrahedron = priv_findFeasiblePointSpectrahedron(SpS);
c = SpS.G * c_spectrahedron + SpS.c;
SpS.center.val = c;

% ------------------------------ END OF CODE ------------------------------
