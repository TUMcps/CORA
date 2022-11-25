function res = testLongDuration_zonotope_plus
% testLongDuration_zonotope_plus - unit test function of plus
%
% Syntax:  
%    res = testLongDuration_zonotope_plus
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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Random Tests ---------------------------------------------------------

tol = 1e-9;
dims = 5:5:100;
testsPerDim = 1000;

% box has to be the same as conversion to interval
for d=1:length(dims)
    for test=1:testsPerDim
        % create a random zonotope
        nrOfGens = randi([10,25],1,1);
        c1 = -1+2*rand(dims(d),1);
        G1 = -1+2*rand(dims(d),nrOfGens);
        Z1 = zonotope(c1,G1);
        c2 = -1+2*rand(dims(d),1);
        G2 = -1+2*rand(dims(d),nrOfGens);
        Z2 = zonotope(c2,G2);
        
        % add zonotopes
        Zres = Z1 + Z2;

        % check center
        if ~all(abs(c1+c2 - center(Zres)) < tol)
            res(d,test) = false;
            continue
        end
        
        % check generators
        Gres = generators(Zres);
        G1inGres = nnz(all(ismember(Gres,G1),1));
        G2inGres = nnz(all(ismember(Gres,G2),1));
        
        if size(Gres,2) ~= G1inGres+G2inGres
            res(d,test) = false;
            continue
        end
        
        % both checks ok
        res(d,test) = true;
    end
end

% combine results
res = all(all(res));

if res
    disp('testLongDuration_zonotope_plus successful');
else
    disp('testLongDuration_zonotope_plus failed');
end

%------------- END OF CODE --------------
