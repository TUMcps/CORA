function options = params2options(params, options)
% params2options - copies all struct fields from param describing the model
%    into the options struct for internal usage
%
% Syntax:
%    options = params2options(params, options)
%
% Inputs:
%    params  - model parameters for object
%    options - algorithm parameters for object
%
% Outputs:
%    options - full options struct for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      23-April-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

paramFields = fieldnames(params);
for i=1:length(paramFields)
    % copy params into options
    % note: if same field already in options, it is overwritten
    options.(paramFields{i}) = params.(paramFields{i});
end


end