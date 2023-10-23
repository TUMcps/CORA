function displayIds(id,varName)
% displayIds - Displays the center and generators of a zonotope
%
% Syntax:
%    displayIds(id,varName)
%
% Inputs:
%    id - id vector
%    varName - name of id variable
%
% Outputs:
%    (to console)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/display

% Authors:       Tobias Ladner
% Written:       31-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

idMin = min(id); idMax = max(id);
disp(varName + ": (" + length(id) + " ids)");

if ~isempty(id) && compareMatrices(id',idMin:idMax,0,"equal",true)
    % short version
    fprintf("    (%g:%g)'\n", idMin, idMax);
else
    % disp all ids
    disp(id);
 end

end

% ------------------------------ END OF CODE ------------------------------
