function res = test_spectraShadow_representsa
% test_spectraShadow_representsa - unit test function of representsa
%
% Syntax:
%    res = test_spectraShadow_representsa
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

% --- origin --------------------------------------------------------------

% fully empty case
SpS = spectraShadow.empty(2);
assert(~representsa(SpS,'origin'));

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
SpS = spectraShadow(polytope(A,b));
assert(~representsa(SpS,'origin'));

% 2D, only origin
A = [1 0; 0 1; -1 0; 0 -1]; b = zeros(4,1);
SpS = spectraShadow(polytope(A,b));
assert(representsa(SpS,'origin'));

% 2D, shifted center
SpS = SpS + [0.01; 0];
assert(~representsa(SpS,'origin'));
% ...add tolerance
tol = 0.02;
assert(representsa(SpS,'origin',tol));

% --- point ---------------------------------------------------------------

% fully empty case
SpS = spectraShadow.empty(2);
assert(~representsa(SpS,'point'));

% 2D, fully empty
A = zeros(0,2); b = zeros(0,0);
SpS = spectraShadow(polytope(A,b));
assert(~representsa(SpS,'point'));

% 2D, only origin
A = [1 0; 0 1; -1 0; 0 -1]; b = zeros(4,1);
SpS = spectraShadow(polytope(A,b));
assert(representsa(SpS,'point'));

% 2D, shifted center
SpS = SpS + [0.01; 0];
assert(representsa(SpS,'point'));
% ...add tolerance
tol = 0.02;
assert(representsa(SpS,'origin',tol));

% --- emptySet ------------------------------------------------------------

% empty constructor
SpS = spectraShadow.empty(2);
assert(representsa(SpS,'emptySet'));

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
SpS = spectraShadow(polytope(A,b));
assert(~representsa(SpS,'emptySet'));

% 1D, empty
A = [1; -1]; b = [1; -3];
SpS = spectraShadow(polytope(A,b));
assert(representsa(SpS,'emptySet'));

% 1D, empty
A = [1; -1; 1; 1; -1; -1]; b = [1; -3; 1; 4; 2; 1];
SpS = spectraShadow(polytope(A,b));
assert(representsa(SpS,'emptySet'));

% 2D, polytope encloses origin
A = [2 1; -2 3; -2 -2; 4 1]; b = ones(4,1);
SpS = spectraShadow(polytope(A,b));
assert(~representsa(SpS,'emptySet'));

% 2D, empty, only inequalities
A = [-1 -1; -1 1; 1 0]; b = [-2; -2; -1];
SpS = spectraShadow(polytope(A,b));
assert(representsa(SpS,'emptySet'));

% 2D, empty only equalities: x1 == 1, x2 == 1, x1+x2 == 1
Ae = [1 0; 0 1; 1 1]; be = [1;1;1];
SpS = spectraShadow(polytope([],[],Ae,be));
assert(representsa(SpS,'emptySet'));

% 2D, unbounded, degenerate (line)
A = [0 1;0 -1]; b = [1;-1];
SpS = spectraShadow(polytope(A,b));
assert(~representsa(SpS,'emptySet'));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
