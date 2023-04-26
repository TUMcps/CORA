function [paramsList,optionsList] = outputParamsOptionsLists
% outputParamsOptionsLists - rewrites the lists for params / options
%    (previously global variables) as a local output
%
% Syntax:
%    [paramsList,optionsList] = outputParamsOptionsLists
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      08-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% read global variables
global fullParamsList;
global fullOptionsList;

% write global variables to output arguments
paramsList = fullParamsList;
optionsList = fullOptionsList;

% delete global variables
clear fullParamsList;
clear fullOptionsList;

end

%------------- END OF CODE --------------
