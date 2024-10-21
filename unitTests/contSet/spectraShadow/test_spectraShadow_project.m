function res = test_spectraShadow_project
% test_spectraShadow_project - unit test function for projection of
%    spectrahedral shadows
%
% Syntax:
%    res = test_spectraShadow_project()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       07-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D, bounded
A = [-1 0; 2 4; 1 -2];
b = [-1; 14; -1];
SpS = spectraShadow(polytope(A,b));
% project to dimension 1, transform to interval, and compare to true result
I_SpS = interval(project(SpS,1));
V = vertices(I_SpS);
V_ = [1 3];
assert(compareMatrices(V,V_,1e-5));
% same with dimension 2
I_SpS = interval(project(SpS,2));
V = vertices(I_SpS);
V_ = [1 3];
assert(compareMatrices(V,V_,1e-5));

% 2D, unbounded
SpS = spectraShadow(polytope([1 0; -1 0; 0 1], [1;1;1]));
% project to dimension 1
SpS_proj = project(SpS,1);
% check if it is still unbounded
assert(isBounded(SpS_proj));

% 2D, degenerate
SpS = spectraShadow(polytope([1 0; -1 0; 0 1; 0 -1], [1;1;1;-1]));
SpS_proj = project(SpS,1);
% spectraShadow should remain unchanged but now full dimensional dimension
assert(isFullDim(SpS_proj));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
