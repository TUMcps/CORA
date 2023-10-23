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
%    Z = zonotope(rand(2,50));
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
    fprintf("    Remainder of " + varName + " (entries " + num2str(maxGens+1) + ...
        "-" + nrOfGens + ") not shown. Check workspace.\n\n");
end

% ------------------------------ END OF CODE ------------------------------
