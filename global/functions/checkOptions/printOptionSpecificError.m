function msg = printOptionSpecificError(obj, property, errmsg)
% printOptionSpecificError - used in checkOptionsReach, checkOptionsSimulate
%    to report a specific error message concerning property
%
% Syntax:  
%    printOptionMissing(obj, property)
%
% Inputs:
%    obj      - system object
%    property - name of missing option
%    errmsg   - error message
%
% Outputs:
%    msg - error message
%
% Example: 
%
% 
% Author:       Mark Wetzlinger
% Written:      03-May-2020
% Last update:  05-May-2020 (output msg)
% Last revision:---


%------------- BEGIN CODE --------------

msg = sprintf(...
    'Error in options check for %s object:\n  options.%s: %s',...
    class(obj),property,errmsg);

end

%------------- END CODE ----------------

