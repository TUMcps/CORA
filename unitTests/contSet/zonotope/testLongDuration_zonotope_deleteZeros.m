function res = testLongDuration_zonotope_deleteZeros
% test_deleteZeros - unit test function of deleteZeros
%    this encompasses checking the function nonzeroFilter
%
% Syntax:  
%    res = testLongDuration_zonotope_deleteZeros
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
% Last update:  09-August-2020 (MW, enhance randomness)
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Random Tests ---------------------------------------------------------

tol = 1e-9;
dims = 5:5:25;
testsPerDim = 1000;

% compare randomly generated zonotopes
for d=1:length(dims)
    for test=1:testsPerDim
        % create random zonotope
        nrOfGens = randi([10,100],1,1);
        
        c = zeros(dims(d),1);
        Gnozeros = -1+2*rand(dims(d),nrOfGens);
        Znozeros = zonotope(c,Gnozeros);
        
        Zdel1 = deleteZeros(Znozeros);
        
        % since no zero generators, results has to be the same as before
        res_rand1 = all(all(abs(Znozeros.Z - Zdel1.Z) < tol));
        
        % insert zero generators in matrix at random position
        Gzeros = [Gnozeros,zeros(dims(d),randi([5,25],1,1))];
        Gzeros = Gzeros(:,randperm(size(Gzeros,2)));
        Zzeros = zonotope(c,Gzeros);
        
        % delete zero generators
        Zdel2 = deleteZeros(Zzeros);
        
        % result has to be the same as original zonotope
        res_rand2 = isequal(Znozeros,Zdel2);

        % add results
        res_rand(d,test) = res_rand1 && res_rand2;
        
    end
end


% add results
res = all(all(res_rand));

if res
    disp('test_deleteZeros successful');
else
    disp('test_deleteZeros failed');
end

%------------- END OF CODE --------------
