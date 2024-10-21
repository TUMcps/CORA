function res = testFlaky_spectraShadow_randPoint
% testFlaky_spectraShadow_randPoint - unit test function of randPoint
%
% Syntax:
%    res = testFlaky_spectraShadow_randPoint
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
% Last update:   04-October-2024 (TL, made flaky)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-5;

% empty set
SpS = spectraShadow.empty(2);
p = randPoint(SpS);
assert(isnumeric(p) && isempty(p) && all(size(p) == [2,0]));

% 1D, bounded
A = [1;-1]; b = [4;-2];
SpS = spectraShadow(polytope(A,b));
p = randPoint(SpS);
p(:,end+1) = randPoint(SpS,1,'extreme');
assert(all(contains(SpS,p,'exact',tol)));

% 1D, single point
Ae = 2; be = 4;
SpS = spectraShadow(polytope([],[],Ae,be));
p = randPoint(SpS);
assert(withinTol(p,2));

% 2D, bounded
A = [2 1; -1 2; -2 -3; 3 -1]; b = ones(4,1);
SpS = spectraShadow(polytope(A,b));
% sample random point with different syntaxes
p = randPoint(SpS);
p(:,end+1) = randPoint(SpS,1);
p(:,end+1) = randPoint(SpS,1,'standard');
p(:,end+1) = randPoint(SpS,1,'uniform');
p(:,end+1) = randPoint(SpS,1,'uniform:hitAndRun');
p(:,end+1) = randPoint(SpS,1,'uniform:ballWalk');
p(:,end+1) = randPoint(SpS,1,'uniform:billiardWalk');
p(:,end+1) = randPoint(SpS,1,'extreme');
p(:,end+1:end+5) = randPoint(SpS,5,'standard');
p(:,end+1:end+5) = randPoint(SpS,5,'uniform');
p(:,end+1:end+5) = randPoint(SpS,5,'uniform:hitAndRun');
p(:,end+1:end+5) = randPoint(SpS,5,'uniform:ballWalk');
p(:,end+1:end+5) = randPoint(SpS,5,'uniform:billiardWalk');
p(:,end+1:end+5) = randPoint(SpS,5,'extreme');
% all have to be contained in P
assert(all(contains(SpS,p,'exact',tol)));

% 2D, bounded, degenerate
A = [1 0; -1 0; 0 1; 0 -1]; b = [1;1;1;-1];
SpS = spectraShadow(polytope(A,b));
p = randPoint(SpS);
assert(contains(SpS,p,'exact',tol));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
