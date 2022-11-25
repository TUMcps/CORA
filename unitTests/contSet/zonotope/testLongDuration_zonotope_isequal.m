function res = testLongDuration_zonotope_isequal
% testLongDuration_zonotope_isequal - unit test function of isequal
%
% Syntax:  
%    res = testLongDuration_zonotope_isequal
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  21-April-2020
%               09-August-2020 (enhance randomness)
% Last revision:---

%------------- BEGIN CODE --------------


% 1. Random Tests ---------------------------------------------------------

dims = 2:8;
testsPerDim = 1000;

% compare randomly generated zonotopes
for d=1:length(dims)
    for test=1:testsPerDim
        % create two random zonotopes
        nrOfGens = randi([10,25],1,1);
        Z1 = zonotope(-1+2*rand(dims(d),nrOfGens+1));
        Z2 = zonotope(-1+2*rand(dims(d),nrOfGens+1));

        % check all combinations
        res_randcomp(d,test) = isequal(Z1,Z1) && isequal(Z2,Z2) && ~isequal(Z1,Z2);
        
        % check different order of generators
        c3 = zeros(dims(d),1);
        G3 = -1+2*rand(dims(d),nrOfGens);
        G3difforder = G3(:,randperm(nrOfGens));
        
        % init zonotopes
        Z3 = zonotope(c3,G3);
        Z3difforder = zonotope(c3,G3difforder);
        
        % compare zonotopes
        res_randorder(d,test) = isequal(Z3,Z3difforder);
    end
end


% add results
res = all(all(res_randcomp)) && all(all(res_randorder));

if res
    disp('test_isequal successful');
else
    disp('test_isequal failed');
end

%------------- END OF CODE --------------
