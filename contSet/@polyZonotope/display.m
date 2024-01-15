function display(pZ)
% display - displays the properties of a polyZonotope object (center,
%    dependent generator matrix, exponent matrix, independent generator
%    matrix, identifiers) on the command window
%
% Syntax:
%    display(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    ---
%
% Example: 
%    pZ = polyZonotope([2;1],[1 0; -2 1],[1; 0],[0 2; 1 0])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       02-May-2020
% Last update:   09-June-2020 (show values)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases
if representsa(pZ,'emptySet')
    dispEmptySet(pZ,inputname(1));
    return
elseif representsa(pZ,'fullspace')
    dispRn(pZ,inputname(1));
    return
end


fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display center
disp('c: ');
disp(center(pZ));

% display dependent generators
displayGenerators(pZ.G,DISPLAYDIM_MAX,'G');

% display independent generators
displayGenerators(pZ.GI,DISPLAYDIM_MAX,'GI');

% display exponential matrix
displayGenerators(pZ.E,DISPLAYDIM_MAX,'E');

% display id
displayIds(pZ.id,'id');

% ------------------------------ END OF CODE ------------------------------
