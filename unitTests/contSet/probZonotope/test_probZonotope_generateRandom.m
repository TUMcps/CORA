function res = test_probZonotope_generateRandom
% test_probZonotope_generateRandom - unit test function of generateRandom
%
% Syntax:
%    res = test_probZonotope_generateRandom
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
% Written:       02-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% test options
Z = probZonotope.generateRandom();
resvec(end+1) = true;

Z = probZonotope.generateRandom('Dimension',3);
resvec(end+1) = true;

Z = probZonotope.generateRandom('Center',ones(2,1));
resvec(end+1) = true;

Z = probZonotope.generateRandom('Dimension',4,'NrGenerators',10);
resvec(end+1) = true;

Z = probZonotope.generateRandom('Distribution','gamma');
resvec(end+1) = true;

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
