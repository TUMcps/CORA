function list = rmiffield(list,field)
% rmiffield - removes a field from a struct if exists
%
% Syntax:
%    str = rmiffield(list,field)
%
% Inputs:
%    list - struct
%    field - field name
%
% Outputs:
%    list - struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       17-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isfield(list,field)
    list = rmfield(list,field);
end

% ------------------------------ END OF CODE ------------------------------
