function [res,msg] = c_partition(val,sys,options)
% c_partition - costum validation function for options.partition
%
% Syntax:
%    [res,msg] = c_partition(val,sys,list)
%
% Inputs:
%    val - value for given param / option
%    sys - linearSys object
%    options - algorithm parameters
%
% Outputs:
%    res - logical whether validation was successful
%    msg - error message if validation failed
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Authors:       Mark Wetzlinger
% Written:       04-March-2019
% Last update:   21-April-2020 (split in lin/poly)
%                03-May-2020 (rewriting of error msgs using class(obj))
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume check ok
res = true;
msg = '';

% check if values in partition are ok
if size(options.partition,2) ~= 2
    % columns have to be [start, end]
    res = false;
    msg = 'has to have 2 columns';
    return;
else
    % check if partition does not overlap, covers all indices
    rows = size(options.partition,1);
    if options.partition(1,1) ~= 1 || ...
            options.partition(rows,2) ~= sys.dim
        res = false;
        msg = '''s indices must not overlap';
        return;
    end
    for i = 2:rows
        if options.partition(i-1,2)+1 ~= options.partition(i,1)
            res = false;
            msg = '''s indices must not leave any dimension out';
            return;
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
