function display(matZ)
% display - Displays the center and generators of a matrix zonotope
%
% Syntax:
%    display(matZ)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    ---
%
% Example: 
%    matZ = matZonotope()
%    C = [0 0; 0 0];
%    G{1} = [1 3; -1 2]; G{2} = [2 0; 1 -1];
%    matZ = matZonotope(C,G)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-June-2010
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(matZ)
    
    dispEmptyObj(matZ,inputname(1));

else

    %display dimension, generators
    disp('dimension: ');
    disp(dim(matZ));
    disp('nr of generators: ');
    disp(matZ.gens);
    %display center
    disp('center: ');
    disp(matZ.center);
    
    %display generators
    disp('generators: ');
    for i=1:length(matZ.generator)
        disp(matZ.generator{i}); 
        disp('---------------'); 
    end

end

% ------------------------------ END OF CODE ------------------------------
