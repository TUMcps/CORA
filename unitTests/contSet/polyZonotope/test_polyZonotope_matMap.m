function res = test_polyZonotope_matMap
% test_polyZonotope_matMap - unit test function of matMap
%
% Syntax:
%    res = test_polyZonotope_matMap
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

% Authors:       Tobias Ladner
% Written:       07-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% settings
nrTest = 10; nrSamples = 100;

% run tests
for i = 1:nrTest
    n = 2; k = 3; m = 2;

    % gather intervals
    I1 = interval.generateRandom('Dimension',n*k);
    I2 = interval.generateRandom('Dimension',k*m);

    % convert to polyZonotopes
    pZ1 = polyZonotope(I1);
    pZ2 = polyZonotope(I2);
    % make them independent
    pZ2 = pZ2.replaceId(max(pZ1.id)+pZ2.id);

    % compute matrix map
    pZ = matMap(pZ1,pZ2,n,k,m);

    % compute samples  
    S1 = reshape(I1.randPoint(nrSamples),[n,k,nrSamples]);
    S2 = reshape(I2.randPoint(nrSamples),[k,m,nrSamples]);
    S = reshape(pagemtimes(S1,S2),[n*m,nrSamples]);

    % check containment
    I = interval(pZ);
    assertLoop(all(contains(I,S)),i)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
