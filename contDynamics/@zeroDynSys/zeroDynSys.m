classdef zeroDynSys < contDynamics
% zeroDynSys class (zeroDynSys: zeor dynamics system; dummy system)
%
% Syntax:  
%    object constructor: Obj = linearSys(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    name - name of system
%    dim - dimension of system
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
% Written:      02-November-2007  
% Last update:  20-March-2008
%               06-August-2016 (changed to new OO format)
% Last revision:---

%------------- BEGIN CODE --------------

methods
    %class constructor
    function obj = zeroDynSys(varargin)
        obj@contDynamics(varargin{1},ones(varargin{2},1),1,1); %instantiate parent class
    end
end
end



%------------- END OF CODE --------------