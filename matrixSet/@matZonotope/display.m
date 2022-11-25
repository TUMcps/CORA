function display(matZ)
% display - Displays the center and generators of a matrix zonotope
%
% Syntax:  
%    display(matZ)
%
% Inputs:
%    matZ - matrix zonotope object
%
% Outputs:
%    ---
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      18-June-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%display dimension, generators
disp('dimension: ');
disp(matZ.dim);
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

%------------- END OF CODE --------------