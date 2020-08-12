function msg = printOptionOutOfRange(obj, property, strct)
% printOptionOutOfRange - used in checkOptionsReach, checkOptionsSimulate
%    to tell if a specified option is out of range
%
% Syntax:  
%    printOptionOutOfRange(obj, property)
%
% Inputs:
%    obj      - system object
%    property - name of wrongly specified option
%    strct    - which struct: params / options
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
%               08-May-2020 (add valid range)
% Last revision:---


%------------- BEGIN CODE --------------

validrange = getValidRange(property);

msg = sprintf(...
    'Error in params/options check for %s object:\n  %s.%s is out of range.\n  valid range: %s.',...
    class(obj),strct,property,validrange);

end

%------------- END CODE ----------------
