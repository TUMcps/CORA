function res = testLong_zonotope_box
% testLong_zonotope_box - unit test function of box
%
% Syntax:
%    res = testLong_zonotope_box
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
% Written:       26-August-2019
% Last update:   09-August-2020 (enhance randomness of test)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 100;

% box has to be the same as conversion to interval
for i=1:nrTests

    % random dimension
    n = randi(20);

    % create a random zonotope
    nrOfGens = 5*n;
    Z = zonotope(-1+2*rand(n,nrOfGens+1));

    % compute axis-aligned box
    Zbox = box(Z);
    cbox = center(Zbox);
    Gbox = generators(Zbox);

    % convert to interval and back to zonotope
    Zint = zonotope(interval(Z));
    cint = center(Zint);
    Gint = generators(Zint);

    % check if axis-aligned box same as interval
    if ~compareMatrices(cbox,cint,1e-14) || ~compareMatrices(Gbox,Gint,1e-14)
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
