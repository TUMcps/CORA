function res = testLong_zonotope_center
% testLong_zonotope_center - unit test function of center
%
% Syntax:
%    res = testLong_zonotope_center
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
% Written:       09-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 1000;

% box has to be the same as conversion to interval
for i=1:nrTests

    % random dimension
    n = randi(20);

    % create a random zonotope
    nrGens = randi([2*n,5*n]);
    c = rand(n,1);
    Z = zonotope.generateRandom('Dimension',n,'Center',c,...
        'NrGenerators',nrGens);

    % compute center
    Zcenter = center(Z);

    % check if centers are the same
    if ~compareMatrices(c,Zcenter)
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
