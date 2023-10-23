classdef setproperty < handle
% setproperty - object constructor for set properties
%
% Syntax:
%    prop = setproperty()
%    prop = setproperty(val)
%
% Inputs:
%    val - value for property
%
% Outputs:
%    obj - generated setproperty object
%
% Example: 
%    val = true;
%    empty = setproperty(val);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       12-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
properties
    % value for property (name of object)
    val
end

methods
    function obj = setproperty(varargin)
        if nargin == 0
            obj.val = [];
        elseif nargin == 1
            obj.val = varargin{1};
        else
            throw(CORAerror('CORA:tooManyInputArgs',2));
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
