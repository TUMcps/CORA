function res = test_zonotope_project
% test_zonotope_project - unit test function of project
%
% Syntax:  
%    res = test_zonotope_project
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

% 1. Analytical Test ------------------------------------------------------

% create zonotope
Z = zonotope([-4, -3, -2, -1; 1, 2, 3, 4; 5, 5, 5, 5]);

% obtain result
Zres = project(Z,[1 3]);

% obtain zonotope matrix
Zmat = Zres.Z;

% true result
true_mat = [-4, -3, -2, -1; ...
            5, 5, 5, 5];

% check result
res_val = all(all(Zmat == true_mat));

% 2. Random Tests ---------------------------------------------------------

dims = 5:5:100;
testsPerDim = 1000;

% box has to be the same as conversion to interval
for d=1:length(dims)
    for test=1:testsPerDim
        % create a random zonotope
        nrOfGens = randi([10,25],1,1);
        c = -1+2*rand(dims(d),1);
        G = -1+2*rand(dims(d),nrOfGens);
        Z = zonotope(c,G);

        % choose random subspace
        projDims = randi([1,dims(d)],1,2);

        % project original center and generator matrix
        cproj = c(projDims);
        Gproj = G(projDims,:);
        
        % project zonotope
        Zproj = project(Z,projDims);
        
        % check projections return the same result
        res_rand(d,test) = all(~any(abs(Zproj.Z - [cproj,Gproj])));
    end
end


% add results
res = res_val && all(all(res_rand));

if res
    disp('test_zonotope_project successful');
else
    disp('test_zonotope_project failed');
end

%------------- END OF CODE --------------
