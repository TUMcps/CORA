function res = test_zonotope_box
% test_zonotope_box - unit test function of box
%
% Syntax:  
%    res = test_zonotope_box
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
% Written:      26-August-2019
% Last update:  09-August-2020 (enhance randomness of test)
% Last revision:---

%------------- BEGIN CODE --------------

% set tolerance
tol = 1e-9;

% 1. Analytical Test ------------------------------------------------------

% init zonotope
Z = zonotope([1;0],[2 -1; 4 1]);

% compute axis-aligned box
Zbox = box(Z);

% convert to interval and back to zonotope
Ztrue = zonotope([1;0],[3 0; 0 5]);

% check if axis-aligned box same as interval
res_val = all(all(abs(Zbox.Z - Ztrue.Z) < tol));

% 2. Random Tests ---------------------------------------------------------

dims = 2:8;
testsPerDim = 1000;

% box has to be the same as conversion to interval
for d=1:length(dims)
    for test=1:testsPerDim
        % create a random zonotope
        nrOfGens = 10;
        Z = zonotope(-1+2*rand(dims(d),nrOfGens+1));

        % compute axis-aligned box
        Zbox = box(Z);

        % convert to interval and back to zonotope
        ZInt = zonotope(interval(Z));

        % check if axis-aligned box same as interval
        res_rand(d,test) = all(all(abs(Zbox.Z - ZInt.Z) < tol));
    end
end

% add results
res = res_val && all(all(res_rand));

if res
    disp('test_box successful');
else
    disp('test_box failed');
end

%------------- END OF CODE --------------