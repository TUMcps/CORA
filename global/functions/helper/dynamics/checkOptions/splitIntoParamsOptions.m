function [params,options_] = splitIntoParamsOptions(options)
% splitIntoParamsOptions - split an options struct (where the parameters
%    and options have been merged) back into a params struct and an options
%    struct; this is necessary in hybrid systems reachability since the
%    options validation merged these, but location/reach has to call
%    contDynamics/reach with separate params/options
%
% Syntax:
%    [params,options_] = splitIntoParamsOptions(options)
%
% Inputs:
%    options - settings for hybrid reachability analysis
%
% Outputs:
%    params - parameters for continuous system reachability analysis
%    options_ - settings for continuous system reachability analysis
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/reach

% Authors:       Mark Wetzlinger
% Written:       21-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init structs
params = struct();
options_ = struct();

% read out all fields of options
allFields = fields(options);

% loop over all fields of options
for i=1:length(allFields)
    % read out name of option
    fieldname = allFields{i};

    % special handling if field is a struct
    if ~isstruct(options.(fieldname))
        % regular option
        if isparam(fieldname)
            params.(fieldname) = options.(fieldname);
        else
            options_.(fieldname) = options.(fieldname);
        end
    else
        % loop over all fields of that option
        subAllFields = fields(options.(fieldname));
        for j=1:length(subAllFields)
            % read out name of option
            subfieldname = subAllFields{j};

            % all nested ones are options
            options_.(fieldname).(subfieldname) = options.(fieldname).(subfieldname);
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
