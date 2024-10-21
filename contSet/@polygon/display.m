function display(pgon)
% display - displays the polygon in the command window
%
% Syntax:
%    display(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases
if representsa(pgon,'emptySet')
    dispEmptySet(pgon,inputname(1));
    return
elseif representsa(pgon,'fullspace')
    dispRn(pgon,inputname(1));
    return
end

% display input variable
fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

%display dimension
display@contSet(pgon);

%display vertices
fprintf('- vertices: %i \n',size(vertices_(pgon),2));
fprintf('- number of regions: %i \n',pgon.set.NumRegions);
fprintf('- number of holes: %i \n',pgon.set.NumHoles);
fprintf(newline);

end

% ------------------------------ END OF CODE ------------------------------
