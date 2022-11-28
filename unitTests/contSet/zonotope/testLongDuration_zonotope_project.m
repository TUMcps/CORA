function res = testLongDuration_zonotope_project
% testLongDuration_zonotope_project - unit test function of project
%
% Syntax:  
%    res = testLongDuration_zonotope_project
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

dims = 5:5:100;
testsPerDim = 1000;

for i=1:length(dims)
    for test=1:testsPerDim
        % create a random zonotope
        nrOfGens = randi([10,25],1,1);
        c = -1+2*rand(dims(i),1);
        G = -1+2*rand(dims(i),nrOfGens);
        Z = zonotope(c,G);

        % choose random subspace
        projDims = randi([1,dims(i)],1,2);
        if rand < 0.5
            % choose logical indexing
            temp = projDims;
            projDims = false(dims(i),1);
            projDims(temp) = true;
        end

        % project original center and generator matrix
        cproj = c(projDims);
        Gproj = G(projDims,:);
        
        % project zonotope
        Zproj = project(Z,projDims);
        
        % check projections return the same result
        res_rand(i,test) = all(~any(abs(Zproj.Z - [cproj,Gproj])));
    end
end


% add results
res = all(all(res_rand));

if ~res
    path = pathFailedTests(mfilename());
    save(path,'nrOfGens','c','G','projDims');
end

%------------- END OF CODE --------------
