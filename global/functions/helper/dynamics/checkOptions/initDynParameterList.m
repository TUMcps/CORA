function [paramsList,optionsList] = initDynParameterList()
% initDynParameterList - inits the struct for dynamic parameter blueprint
%
% Syntax:
%    [paramsList,optionsList] = initDynParameterList()
%
% Inputs:
%    -
%
% Outputs:
%    paramsList - struct
%    optionsList - struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: config files, add2list

% Authors:       Tobias Ladner
% Written:       06-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

paramsList = struct('name',{},'status',{},'condfun',{});
optionsList = struct('name',{},'status',{},'condfun',{});

% ------------------------------ END OF CODE ------------------------------
