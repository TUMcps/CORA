function printDefaultValue(obj, property, defValue)
% printDefaultValue - used in checkOptionsReach, checkOptionsSimulate
%    to tell if a specified option is equivalent to its default setting
%    current status: no warning, function kept for possible future use
%
% Syntax:  
%    printDefaultValue(obj, property, defValue)
%
% Inputs:
%    obj      - system object
%    property - name of option
%    defValue - default value of option
%
% Outputs:
%    -
%
% Example: 
%
% 
% Author:       Mark Wetzlinger
% Written:      08-Aug-2019
% Last update:  03-May-2020
%               19-May-2020 (no warning for def values)
% Last revision:---


%------------- BEGIN CODE --------------

% warning("Options check for " + class(obj) + ": default value for option " + property + " is " + defValue + ...
%        ", therefore setting it explicitly to " + defValue + " is redundant");

end

%------------- END CODE ----------------

