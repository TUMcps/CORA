function res = testLongDuration_zonotope_center
% testLongDuration_zonotope_center - unit test function of center
%
% Syntax:  
%    res = testLongDuration_zonotope_center
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

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020 (MW, enhance randomness of test)
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% number of tests
nrTests = 1000;

% box has to be the same as conversion to interval
for i=1:nrTests

    % random dimension
    n = randi(20);

    % create a random zonotope
    nrOfGens = randi([2*n,5*n]);
    c = rand(n,1);
    Z = zonotope(c,-1+2*rand(n,nrOfGens));

    % compute center
    Zcenter = center(Z);

    % check if centers are the same
    if ~compareMatrices(c,Zcenter)
        res = false; return
    end
end

%------------- END OF CODE --------------
