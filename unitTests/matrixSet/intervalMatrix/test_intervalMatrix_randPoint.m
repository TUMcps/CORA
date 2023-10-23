function res = test_intervalMatrix_randPoint
% test_intervalMatrix_randPoint - unit test function for random sampling
% 
% Syntax:
%    res = test_intervalMatrix_randPoint
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

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate interval matrix
intMat = intervalMatrix([2 3; 1 2],[1 0; 1 1]);

% sample random matrix
M = randPoint(intMat);

% check containment
res = contains(intMat,M{1});

% multiple samples
M = randPoint(intMat,5);
for i=1:length(M)
    if ~contains(intMat,M{i})
        res = false; break
    end
end

% extreme sampling
M = randPoint(intMat,1,'extreme');
res = res & contains(intMat,M{1},'exact',1e-8);

% ------------------------------ END OF CODE ------------------------------
