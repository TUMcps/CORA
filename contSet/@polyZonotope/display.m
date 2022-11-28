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

% Author:        Mark Wetzlinger
% Written:       02-May-2020
% Last update:   09-June-2020 (show values)
% Last revision: ---

%------------- BEGIN CODE --------------

if isempty(pZ)
    
    dispEmptyObj(pZ,inputname(1));
    
else
    
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
    
    % display center
    disp('c: ');
    disp(center(pZ));
    
    % display dependent generators
    G = pZ.G;
    displayGenerators(G,DISPLAYDIM_MAX,'G');
    
    % display independent generators
    Grest = pZ.Grest;
    displayGenerators(Grest,DISPLAYDIM_MAX,'Grest');
    
    % display exponential matrix
    expMat = pZ.expMat;
    displayGenerators(expMat,DISPLAYDIM_MAX,'expMat');
    
    % display id
    disp('id:');
    disp(pZ.id);
    
end

%------------- END OF CODE --------------
