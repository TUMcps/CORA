function display(cZ)
% display - Displays the properties of a conZonotope object (center,
%    generators, and constraints for the factors) on the command window
%
% Syntax:
%    display(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    ---
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       20-December-2017
% Last update:   01-May-2020 (MW, added empty case)
%                09-June-2020 (MW, restrict number of shown generators)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases
if representsa(cZ,'emptySet')
    dispEmptySet(cZ,inputname(1));
    return
elseif representsa(cZ,'fullspace')
    dispRn(cZ,inputname(1));
    return
end


fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display center and generators
disp('c: ');
disp(cZ.c);

G = cZ.G;
displayGenerators(G,DISPLAYDIM_MAX,'G');

%display constraint system
if isempty(cZ.A) && isempty(cZ.b)
    disp('Constraint system (Ax = b): no constraints.')
    fprintf(newline);
    
else
    disp('Constraint system (Ax = b):')

    displayGenerators(cZ.A,DISPLAYDIM_MAX,'A');

    disp('b: ');
    disp(cZ.b);
end

end

% ------------------------------ END OF CODE ------------------------------
