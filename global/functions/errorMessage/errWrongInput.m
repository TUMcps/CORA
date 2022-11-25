function [msg,id] = errWrongInput(arg)
% errWrongInput - standardized error message format for wrong input
%                 arguments
%
% Syntax:  
%    [msg,id] = errWrongInput(arg)
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

    msg = sprintf('Wrong value for input argument "%s"!',arg);
    id = "CORA:wrongValue";

%------------- END OF CODE --------------