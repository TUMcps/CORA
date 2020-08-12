function msg = printRedundancies(options, validFields)
% printRedundancies - prints all specified options which are not contained
%   in validFields (validFields contains all mandatory and optional
%   parameters for the system)
%
% Syntax:
%    printRedundancies(options, validFields)
%
% Inputs:
%    options     - options for object
%    validFields - list (cell of strings) of non-redundant options
%
% Outputs:
%    msg - warning message
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
% Written:      05-Mar-2019
% Last update:  ---
% Last revision:24-July-2019
%               08-May-2020 (MW, single warning for all redundancies)

%------------- BEGIN CODE --------------

fields = fieldnames(options);
redundantOptions = fields(~ismember(fields,validFields));

if ~isempty(redundantOptions)
    % some param/option is redundant
    redundantOptionsStr = "";
    for i=1:numel(redundantOptions)-1
        if ~any(strcmp(fields{i},validFields))
            redundantOptionsStr = redundantOptionsStr + redundantOptions{i} + ", ";
        end
    end
    redundantOptionsStr = redundantOptionsStr + redundantOptions{end}; % no ,
    
    msg = sprintf(...
        'The following params/options have been set, but are redundant:\n  %s',...
        redundantOptionsStr);

else
    % no redundancies
    msg = '';

end


end

%------------- END OF CODE --------------
