function res = testLong_spectraShadow_center
% testLong_spectraShadow_center - unit test function of center
%
% Syntax:
%    res = testLong_spectraShadow_center
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
% Written:       15-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

runs = 5;

% Test the center of polytopes as spectrahedra
for i=1:runs
    dimension = randi(8);
    P = polytope.generateRandom('Dimension', dimension);
    
    SpS = spectraShadow(P);
    c = center(SpS);
    
    assert(contains(SpS,c,'exact',1e-5))
end

% Do the same with ellipsoids
for i=1:runs
    dimension = randi(8);
    E = ellipsoid.generateRandom('Dimension', dimension);
    
    SpS = spectraShadow(E);
    c = center(SpS);
    
    assert(contains(SpS,c,'exact',1e-5))
end

% ------------------------------ END OF CODE ------------------------------
