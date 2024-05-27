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
%    G = []; G(:,:,1) = [1 3; -1 2]; G(:,:,2) = [2 0; 1 -1];
%    matZ = matZonotope(C,G)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-June-2010
% Last update:   25-April-2024 (TL, harmonized display with contSet classes)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(matZ)
    dispEmptyObj(matZ,inputname(1));
    return
end

% display input variable
fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display dimension
fprintf('%s:\n', class(matZ))
disp('dimension: ');
dims = dim(matZ);
disp(dims);

% compute number of dims to display
dispdims = min(dims,DISPLAYDIM_MAX);
dispgens = min(matZ.numgens,DISPLAYDIM_MAX);

% display center
disp('center: ');
disp(matZ.C(1:dispdims(1),1:dispdims(2)));

% display generators
fprintf('generators: (%i generators)\n', numgens(matZ));
disp(matZ.G(1:dispdims(1),1:dispdims(2),1:dispgens))

addEmptyLine = false;
if any(dispdims < dims)
    fprintf("    Remainder of dimensions (>%i) not shown. Check workspace.\n", DISPLAYDIM_MAX);
    addEmptyLine = true;
end
if any(dispgens < matZ.numgens)
    fprintf("    Remainder of generators (>%i) not shown. Check workspace.\n", DISPLAYDIM_MAX);
    addEmptyLine = true;
end
if addEmptyLine
    fprintf('\n')
end

if isempty(matZ.G)
    fprintf('\n')
end

end

% ------------------------------ END OF CODE ------------------------------
