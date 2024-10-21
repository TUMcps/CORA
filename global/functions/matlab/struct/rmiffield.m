function list = rmiffield(list,field)
% rmiffield - removes a field from a struct if exists
%
% Syntax:
%    list = rmiffield(list,field)
%
% Inputs:
%    list - struct
%    field - field name (cell of char arrays, char array)
%
% Outputs:
%    list - struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner, Mark Wetzlinger
% Written:       17-October-2023
% Last update:   30-August-2024 (MW, support cell-array of fields)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if iscell(field)
    % cell array of fields
    for i=1:length(field)
        if isfield(list,field{i})
            list = rmfield(list,field{i});
        end 
    end
else
    % single field
    if isfield(list,field)
        list = rmfield(list,field);
    end
end

% ------------------------------ END OF CODE ------------------------------
