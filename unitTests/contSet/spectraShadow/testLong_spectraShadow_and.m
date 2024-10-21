function res = testLong_spectraShadow_and
% testLong_spectraShadow_and - unit test function of and
%
% Syntax:
%    res = testLong_spectraShadow_and
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

res = true;

runs = 2;
samples = 5;

% Test the intersection of polytopes as spectrahedra
for i=1:runs
    dimension = randi(8);
    P1 = polytope.generateRandom('Dimension', dimension);
    P2 = polytope.generateRandom('Dimension', dimension);
    % Make sure that some sort of intersection is happening
    P1 = P1 - center(P1);
    P2 = P2 - center(P2);
    
    P = P1 & P2;
    SpS = spectraShadow(P1) & spectraShadow(P2);
    
    p_P = randPoint(P,samples);
    p_SpS = randPoint(SpS,samples);
    
    for j=1:samples
        assertLoop(P.contains(p_SpS(:,j),'exact',1e-6),i,j);
        assertLoop(SpS.contains(p_P(:,j),'exact',1e-6),i,j);
    end
end

% Do the same with ellipsoids, but only in one direction (the exact
% intersection of two ellipsoids can not be represented as an ellipsoid)
for i=1:runs
    dimension = randi(8);
    E1 = polytope.generateRandom('Dimension', dimension);
    E2 = polytope.generateRandom('Dimension', dimension);
    % Make sure that some sort of intersection is happening
    E1 = E1 - center(E1);
    E2 = E2 - center(E2);
    
    E = E1 & E2;
    SpS = spectraShadow(E1) & spectraShadow(E2);
    
    p_SpS = randPoint(SpS,samples);
    
    for j=1:samples
        assertLoop(E.contains(p_SpS(:,j),'exact',1e-6),i,j)
    end
end

% Double check that empty-set computation works
A = [-1 0 0];
SpS = spectraShadow(A);
representsa(SpS,'emptySet');
SpS = SpS & SpS;
assert(~isempty(SpS.emptySet.val) && (isempty(SpS.emptySet.val) || SpS.emptySet.val))

% ------------------------------ END OF CODE ------------------------------
