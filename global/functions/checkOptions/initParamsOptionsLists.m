function initParamsOptionsLists
% initParamsOptionsLists - initializes the lists for params / options
%    as global variables for easier access in validateOptions
%
% Syntax:
%    initParamsOptionsLists
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

% define as global variables for simpler syntax in config file
global fullParamsList;
global fullOptionsList;

% set to empty struct
fullParamsList = struct();
fullOptionsList = struct();

% introduce fields for params list
fullParamsList.name = {};
fullParamsList.status = {};
fullParamsList.checkfuncs = {};
fullParamsList.errmsgs = {};
fullParamsList.condfunc = {};

% introduce fields for options list
fullOptionsList.name = {};
fullOptionsList.status = {};
fullOptionsList.checkfuncs = {};
fullOptionsList.errmsgs = {};
fullOptionsList.condfunc = {};

end

%------------- END OF CODE --------------