function check_enclose(options, obj)
% check_enclose - checks if options.enclose
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    check_enclose(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    -
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
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

option = 'enclose';
validEnclose = {'box','pca','flow'};
strct = 'options';

% enclose: string with values 'box', 'pca' or 'flow'
% not needed if guardIntersect = pancake / hyperplaneMap
if ~isfield(options,option) && ...
        ~strcmp(options.guardIntersect, 'pancake') && ...
        ~strcmp(options.guardIntersect, 'hyperplaneMap') && ...
        ~strcmp(options.guardIntersect, 'levelSet')
    
    error(printOptionMissing(obj,option,strct));
    
elseif ~strcmp(options.guardIntersect, 'pancake') && ...
       ~strcmp(options.guardIntersect, 'hyperplaneMap') && ...
       ~strcmp(options.guardIntersect, 'levelSet')
    
    % check if it is a cell-array
    if ~iscell(options.enclose)
        error(printOptionSpecificError(obj,option,...
            'enclosure has to be a cell-array.'));
    end
    
    % check if entries are valid entries
    for i=1:length(options.enclose)
        if ~ismember(options.enclose{i},validEnclose)
            error(printOptionOutOfRange(obj,option,strct));
        end
    end
end


end

%------------- END OF CODE --------------

