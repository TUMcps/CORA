function display(Z)
% display - Displays the properties of a zonotope object (center and
%    generator matrix) on the command window
%
% Syntax:
%    display(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    ---
%
% Example: 
%    Z = zonotope([1;0], [1 2 -1; 0 -1 1]);
%    display(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       14-September-2006 
% Last update:   22-March-2007
%                27-August-2019
%                01-May-2020 (MW, add empty case)
%                09-June-2020 (MW, restrict number of shown generators)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases
if representsa(Z,'emptySet')
    dispEmptySet(Z,inputname(1));
    return
elseif representsa(Z,'fullspace')
    dispRn(Z,inputname(1));
    return
end


fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

%display dimension
display@contSet(Z);
fprintf(newline);

%display center
disp('c: ');
disp(Z.c);

%display generators
displayGenerators(Z.G,DISPLAYDIM_MAX,'G');

% ------------------------------ END OF CODE ------------------------------
