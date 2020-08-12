function msg = printOptionMissing(obj, property, strct)
% printOptionMissing - used in checkOptionsReach, checkOptionsSimulate
%  to tell if a specified option is out of range
%
% Syntax:  
%    printOptionMissing(obj, property)
%
% Inputs:
%    obj      - system object
%    property - name of missing option
%    strct    - which struct: 'params'/'options'
%
% Outputs:
%    msg - error message
%
% Example: 
%
% 
% Author:       Mark Wetzlinger
% Written:      16-Feb-2019
% Last update:  03-May-2020
%               05-May-2020 (output msg)
% Last revision:---


%------------- BEGIN CODE --------------

error(sprintf(...
    'Error in params/options check for %s object:\n  %s.%s is missing.',...
    class(obj),strct,property));

end

%------------- END CODE ----------------

