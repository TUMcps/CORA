classdef linProbSys < contDynamics
% linProbSys class (linProbSys: linear probabilistic system)
%
% Syntax:  
%    object constructor: Obj = linearSys(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    name - name of system
%    A - state matrix
%    B - input matrix
%    C - noise matrix
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      06-October-2007 
% Last update:  26-February-2008
%               05-August-2016 (changed to new OO format)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    A = [];
    B = [];
    C = [];
    taylor = [];
end

methods
    %class constructor
    function obj = linProbSys(varargin)
        obj@contDynamics(varargin{1},length(varargin{2}),1,1); %instantiate parent class
        %4 inputs
        if nargin==4
            obj.A = varargin{2};
            obj.B = varargin{3};
            obj.C = varargin{4};
        end
    end
end
end


%------------- END OF CODE --------------