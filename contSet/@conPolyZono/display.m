function display(cPZ)
% display - displays the properties of a conPolyZono object (center,
%    generator matrices, exponent matrix, constraint system) on the 
%    command window
%
% Syntax:
%    display(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    ---
%
% Example: 
%    c = [0;0];
%    G = [1 0 1 -1; 0 1 1 1];
%    E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    EC = [0 1 2; 1 0 0; 0 1 0];
% 
%    cPZ = conPolyZono(c,G,E,A,b,EC)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/display

% Authors:       Niklas Kochdumper
% Written:       19-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases
if representsa(cPZ,'emptySet')
    dispEmptySet(cPZ,inputname(1));
    return
elseif representsa(cPZ,'fullspace')
    dispRn(cPZ,inputname(1));
    return
end


fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display center
disp('c: ');
disp(cPZ.c);

% display generators
displayGenerators(cPZ.G,DISPLAYDIM_MAX,'G');

% display exponential matrix
displayGenerators(cPZ.E,DISPLAYDIM_MAX,'E');

% display constraint offset
disp('b:');
disp(cPZ.b);

% display constraint generators
displayGenerators(cPZ.A,DISPLAYDIM_MAX,'A');

% display constraint exponential matrix
displayGenerators(cPZ.EC,DISPLAYDIM_MAX,'EC');

% display independent generators
displayGenerators(cPZ.GI,DISPLAYDIM_MAX,'GI');

% display id
displayIds(cPZ.id,'id');

% ------------------------------ END OF CODE ------------------------------
