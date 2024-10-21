function options = setDefaultFields(options,defaultFields)
% setDefaultFields - set default fields for a given options struct.
%
% Syntax:
%    options = nnHelper.setDefaultFields(options,defaultFields)
%
% Inputs:
%    options - options struct
%    name - options name
%    defaultFields - cell containing name and default value
%
% Outputs:
%    options - updated options
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: validateNNoptions, validateRLoptions

% Authors:       Lukas Koller
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Set default value of fields if required.
for i=1:size(defaultFields, 1)
    field = defaultFields{i, 1};
    if ~isfield(options, field)
        fieldValue = defaultFields{i, 2};
        if isa(fieldValue, "function_handle")
            fieldValue = fieldValue(options);
        end
        options.(field) = fieldValue;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
