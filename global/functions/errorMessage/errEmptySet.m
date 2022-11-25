function [msg,id] = errEmptySet()
% errEmptySet - standardized error message format if a set is empty
%
% Syntax:  
%    [msg,id] = errEmptySet()
%
% Inputs:
%    ---
%
% Outputs:
%    msg - error message
%    id - error identifier
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: noops

% Author:        Niklas Kochdumper
% Written:       25-January-2021 
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    msg = "Set is empty!";
    id = "CORA:emptySet";

%------------- END OF CODE --------------
