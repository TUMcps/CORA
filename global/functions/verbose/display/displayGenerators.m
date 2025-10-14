function displayGenerators(G,maxGens,varName)
% displayGenerators - Displays the center and generators of a zonotope
%
% Syntax:
%    displayGenerators(G,maxGens,varName)
%
% Inputs:
%    G - generator matrix
%    maxGens - max number of displayed generators
%    varName - name of generator matrix variable
%
% Outputs:
%    (to console)
%
% Example: 
%    Z = zonotope([ 0.444 ; 0.756 ], [ 0.603 0.114 0.849 0.466 0.630 0.580 0.600 0.035 0.408 0.460 0.551 0.701 0.052 0.460 ; 0.783 0.979 0.051 0.326 0.230 0.603 0.448 0.514 0.108 0.451 0.805 0.872 0.220 0.959 ]);
%    G = generators(Z);
%    displayGenerators(G,10,'G');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/display

% Authors:       Mark Wetzlinger
% Written:       09-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%display generators
nrOfGens = size(G,2);
if nrOfGens == 1
    genStr = "generator";
else
    genStr = "generators";
end
disp(varName + ": (" + nrOfGens + " " + genStr + ")");
if nrOfGens <= maxGens
    disp(G);
else
    disp(G(:,1:maxGens));
    fprintf("    Remainder of generators (>%i) of %s not shown. Check workspace.\n\n", DISPLAYDIM_MAX,varName);
end

% ------------------------------ END OF CODE ------------------------------
