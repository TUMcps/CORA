function res = testLong_spectraShadow_enclose
% testLong_spectraShadow_enclose - unit test function of enclose
%
% Syntax:
%    res = testLong_spectraShadow_enclose
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

% Test enclose of polytopes as spectrahedra
for i=1:runs
    dimension = randi(3);
    M = randn(dimension);
    P1 = polytope.generateRandom('Dimension', dimension);
    P2 = polytope.generateRandom('Dimension', dimension);
    
    SpS = enclose(spectraShadow(P1), M, spectraShadow(P2));
    
    p = randPoint(enclose(P1,M,P2),samples);
    
    
    for j=1:samples
        assertLoop(SpS.contains(p(:,j),'exact',1e-6),i,j);
    end
end


% ------------------------------ END OF CODE ------------------------------
