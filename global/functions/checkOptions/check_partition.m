function check_partition(options, obj)
% check_partition - checks if options.partition
%  1) takes an allowed value
%
% Syntax:
%    check_partition(options, obj)
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
% Written:      26-June-2019
% Last update:  03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

if strcmp(options.linAlg,'decomp')
    option = 'partition';
    strct = 'options';
    if ~isfield(options,option)
        error(printOptionMissing(obj,option,'options'));
    else
        % check if values in partition are ok
        if size(options.partition,2) ~= 2
            % columns have to be [start, end]
            error(printOptionOutOfRange(obj,option,strct));
        else
            % check if partition does not overlap, covers all indices
            rows = size(options.partition,1);
            if options.partition(1,1) ~= 1 || ...
                    options.partition(rows,2) ~= length(center(options.R0))
                error(printOptionOutOfRange(obj,option,strct));
            end
            for i = 2:rows
                if options.partition(i-1,2)+1 ~= options.partition(i,1)
                    error(printOptionOutOfRange(obj,option,strct));
                end
            end
            
        end
    end
end

%------------- END OF CODE --------------

