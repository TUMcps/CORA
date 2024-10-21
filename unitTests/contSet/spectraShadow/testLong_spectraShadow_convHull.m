function res = testLong_spectraShadow_convHull
% testLong_spectraShadow_convHull - unit test function of convHull
%
% Syntax:
%    res = testLong_spectraShadow_convHull
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

% Test the convex hull of polytopes as spectrahedra
for i=1:runs
    dimension = randi(8);
    P1 = polytope.generateRandom('Dimension', dimension);
    P2 = polytope.generateRandom('Dimension', dimension);
    
    SpS = convHull(spectraShadow(P1), spectraShadow(P2));
    
    p_P1 = randPoint(P1,samples);
    p_P2 = randPoint(P2,samples);
    
    lambdas = rand([samples 1]);
    
    for j=1:samples
        assertLoop(SpS.contains(p_P1(:,j) * lambdas(j) + p_P2(:,j) * (1-lambdas(j)),'exact',1e-6),i,j);
    end
end

% Test the convex hull of ellipsoids as spectrahedra
for i=1:runs
    dimension = randi(8);
    E1 = ellipsoid.generateRandom('Dimension', dimension);
    E2 = ellipsoid.generateRandom('Dimension', dimension);
    
    SpS = convHull(spectraShadow(E1), spectraShadow(E2));
    
    p_E1 = randPoint(E1,samples);
    p_E2 = randPoint(E2,samples);
    
    lambdas = rand([samples 1]);
    
    for j=1:samples
        assertLoop(SpS.contains(p_E1(:,j) * lambdas(j) + p_E2(:,j) * (1-lambdas(j)),'exact',1e-6),i,j)
    end
end


% ------------------------------ END OF CODE ------------------------------
