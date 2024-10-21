function res = test_spectraShadow_center
% test_spectraShadow_center - unit test function of center
%
% Syntax:
%    res = test_spectraShadow_center
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
% Written:       05-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% Explicit example
A0 = eye(3);
Ai{1} = [0 1 0;1 0 0;0 0 0];
Ai{2} = [0 0 1;0 0 0;1 0 0];

SpS = spectraShadow([A0,Ai{1},Ai{2}]); %- unit ball

c = center(SpS);

assert(SpS.contains(c, 'exact', 1e-5));

% Construct based on explicit ellipsoid
E = ellipsoid(eye(3));
SpS = spectraShadow(E);
assert(res & (norm(center(SpS))==0));

% ------------------------------ END OF CODE ------------------------------
