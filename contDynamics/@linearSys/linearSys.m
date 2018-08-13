classdef linearSys < contDynamics
% linearSys class (linSys: linear system)
%
% Syntax:  
%    object constructor: Obj = linearSys(varargin)
%                        Obj = linearSys(name,A,B)
%                        Obj = linearSys(name,A,B,c)
%                        Obj = linearSys(name,A,B,c,C)
%                        Obj = linearSys(name,A,B,c,C,D)
%                        Obj = linearSys(name,A,B,c,C,D,k)
%    copy constructor: Obj = otherObj
%
% Description:
%    Generates a linear system object according to the following first
%    order differential equations:
%       x' = A x + B u + c
%       y = C x + D u + k
%
% Inputs:
%    name - name of system
%    A - state matrix
%    B - input matrix
%    c - constant input
%    C - output matrix
%    D - throughput matrix
%    k - output offset
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      23-January-2007 
% Last update:  30-April-2007
%               04-August-2016 (changed to new OO format)
%               01-November-2017 (constant input added)
%               20-March-2018 (NK, output equation parameter added)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    A = []; % system matrix
    B = []; % input matrix
    c = []; % constant input
    C = []; % output matrix
    D = []; % throughput matrix
    k = []; % output offset
    taylor = [];
end

methods
    %class constructor
    function obj = linearSys(varargin)
        obj@contDynamics(varargin{1},ones(length(varargin{2}),1),1,1); %instantiate parent class
        %3 or more inputs
        if nargin>=3
            obj.A = varargin{2};
            obj.B = varargin{3};
            if nargin >= 4
                obj.c = varargin{4};
            end
            if nargin >= 5
                obj.C = varargin{5};
            end
            if nargin >= 6
                obj.D = varargin{6};
            end
            if nargin == 7
                obj.k = varargin{7}; 
            end
        end
    end
end
end

%------------- END OF CODE --------------