function res = testLong_spectraShadow_generateESumRep
% testLong_spectraShadow_generateESumRep - unit test function of generateESumRep
%
% Syntax:
%    res = testLong_spectraShadow_generateESumRep
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

% Test if transforming to ESumRep and back yields the same spectrahedral
% shadow
for i=1:runs
    dimension = randi(8);
    SpS = spectraShadow.generateRandom('Dimension',dimension);
    generateESumRep(SpS);
    
    SpS_copy = spectraShadow(SpS.ESumRep.val);
    
    p_SpS = randPoint(SpS,samples);
    p_SpS_copy = randPoint(SpS_copy,samples);
    
    for j=1:samples
        assertLoop(SpS_copy.contains(p_SpS(:,j),'exact',1e-5) || ~SpS.contains(p_SpS_copy(:,j),'exact',1e-5),i,j)
    end
end


% ------------------------------ END OF CODE ------------------------------
