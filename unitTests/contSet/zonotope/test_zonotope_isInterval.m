function res = test_zonotope_isInterval
% test_zonotope_isInterval - unit test function of isInterval
%
% Syntax:  
%    res = test_zonotope_isInterval
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
% Written:      09-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% 1. Analytical Tests -----------------------------------------------------

% create zonotopes
c1 = [0; 0];
G1 = [2 0; 0 1];
Z1 = zonotope(c1,G1);

c2 = [1; 0];
G2 = [2 1; -1 4];
Z2 = zonotope(c2,G2);

% check result
res_val = isInterval(Z1) && ~isInterval(Z2);


% 2. Random Tests ---------------------------------------------------------

dims = 5:5:100;
testsPerDim = 1000;

% compare randomly generated zonotopes
for d=1:length(dims)
    for test=1:testsPerDim
        % create two random zonotopes
        nrOfGens = randi([10,25],1,1);
        
        c = zeros(dims(d),1);
        GnonInt = -1+2*rand(dims(d),nrOfGens);
        ZnonInt = zonotope(c,GnonInt);
        
        % generate random matrix with only one non-negative entry per column
        GInt = zeros(dims(d),nrOfGens);
        % non-negative row for each generator
        idx = randi([1,dims(d)],1,nrOfGens);
        % linear indexing
        idx = idx + [0:dims(d):dims(d)*(nrOfGens-1)];
        % random values
        vals = -1+2*rand(1,nrOfGens);
        % write random values in matrix
        GInt(idx) = vals;
        % init zonotope
        ZInt = zonotope(c,GInt);

        % check zonotopes
        res_randcomp(d,test) = isInterval(ZInt) && ~isInterval(ZnonInt);
        
    end
end


% add results
res = all(res_val) && all(all(res_randcomp));

if res
    disp('test_isInterval successful');
else
    disp('test_isInterval failed');
end

%------------- END OF CODE --------------
