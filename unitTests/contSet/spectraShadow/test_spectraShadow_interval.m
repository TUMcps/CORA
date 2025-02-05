function res = test_spectraShadow_interval
% test_spectraShadow_interval - unit test function of interval
%
% Syntax:
%    res = test_spectraShadow_interval
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
% See also: -

% Authors:       Adrian Kulmburg
% Written:       14-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


tol = 1e-5;

% empty spectrahedron
SpS_empty = spectraShadow.empty(3);
I_empty = interval(SpS_empty);
assert(representsa(I_empty, 'emptySet'));


% 1D, bounded, non-degenerate
SpS = spectraShadow([1 0 1 0;0 1 0 -1]);
I = interval(SpS);
I_true = interval(-1,1);
assert(isequal(I,I_true,tol));

% 1D, empty
SpS = spectraShadow([-1 0]);
I = interval(SpS);
I_true = interval.empty(1);
assert(isequal(I,I_true,tol));

% 1D, unbounded
SpS = spectraShadow([1 0]);
I = interval(SpS);
I_true = interval(-Inf,Inf);
assert(isequal(I,I_true,tol));

% 1D, single point
SpS = spectraShadow([-1 0 1 0;0 1 0 -1]);
I = interval(SpS);
I_true = interval(1,1);
assert(isequal(I,I_true,tol));


% 2D, bounded, non-degenerate
A0 = eye(4);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
I = interval(SpS);
I_true = interval([-1;-1],[1;1]);
assert(isequal(I,I_true,tol));

% 2D, bounded, degenerate
A0 = blkdiag([-1 0;0 1],[-1 0;0 1]);
Ai{1} = blkdiag([1 0;0 -1],zeros(2));
Ai{2} = blkdiag(zeros(2),[1 0;0 -1]);
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
I = interval(SpS);
I_true = interval([1;1],[1;1]);
assert(isequal(I,I_true,tol));

% 2D, unbounded, non-degenerate 
SpS = spectraShadow([1 0 0]);
I = interval(SpS);
I_true = interval([-Inf;-Inf],[Inf;Inf]);
assert(isequal(I,I_true,tol));

% 2D, unbounded, degenerate
A0 = [-1 0;0 1];
Ai{1} = zeros(2);
Ai{2} = [1 0;0 -1];
SpS = spectraShadow([A0 Ai{1} Ai{2}]);
I = interval(SpS);
I_true = interval([-Inf;1],[Inf;1]);
assert(isequal(I,I_true,tol));

% 2D, empty
SpS = spectraShadow([-1 0 0]);
I = interval(SpS);
I_true = interval.empty(2);
assert(isequal(I,I_true,tol));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
