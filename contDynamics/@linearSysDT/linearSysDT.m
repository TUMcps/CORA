classdef linearSysDT < contDynamics
% linearSysDT class (linSysDT: discrete-time linear system)
%
% Syntax:  
%    obj = linearSysDT(A,B,dt)
%    obj = linearSysDT(A,B,c,dt)
%    obj = linearSysDT(A,B,c,C,dt)
%    obj = linearSysDT(A,B,c,C,D,dt)
%    obj = linearSysDT(A,B,c,C,D,k,dt)
%    obj = linearSysDT(name,A,B,dt)
%    obj = linearSysDT(name,A,B,c,dt)
%    obj = linearSysDT(name,A,B,c,C,dt)
%    obj = linearSysDT(name,A,B,c,C,D,dt)
%    obj = linearSysDT(name,A,B,c,C,D,k,dt)
%
% Description:
%    Generates a discrete-time linear system object according to the 
%    following first-order difference equations:
%       x(i+1) = A x(i) + B u(i) + c + w(i)
%       y(i) = C x(i) + D u(i) + k + v(i)
%
% Inputs:
%    name - name of system
%    dt - sampling time
%    A - state matrix
%    B - input matrix
%    c - constant input
%    C - output matrix
%    D - throughput matrix
%    k - output offset
%
% Outputs:
%    obj - Generated Object
%
% Example:
%    A = [-0.4 0.6; 0.6 -0.4];
%    B = [0; 1];
%    C = [1 0];
%
%    dt = 0.4;
%
%    sys = linearSysDT(A,B,[],C,dt)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      20-Mar-2020 
% Last update:  14-Jun-2021 (MA, invoke observe from superclass)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    A = [];     % system matrix
    B = 1;      % input matrix
    c = [];     % constant input
    C = 1;      % output matrix
    D = [];     % throughput matrix
    k = [];     % output offset
    dt = [];    % sampling time
end

methods
    
    % class constructor
    function obj = linearSysDT(varargin)
        
        % check if system name is provided
        name = 'linearSysDT';
        A = []; B = []; c = []; C = 1; D = []; k = [];
        
        if ischar(varargin{1})
           name = varargin{1};
           A = varargin{2};
           B = varargin{3};
           if nargin == 4
               dt = varargin{4};
           elseif nargin == 5
               c = varargin{4};
               dt = varargin{5};
           elseif nargin == 6
               c = varargin{4};
               C = varargin{5};
               dt = varargin{6};
           elseif nargin == 7
               c = varargin{4};
               C = varargin{5};
               D = varargin{6};
               dt = varargin{7};
           elseif nargin == 8
               c = varargin{4};
               C = varargin{5};
               D = varargin{6};
               k = varargin{7};
               dt = varargin{8};
           else
               error('Wrong number of input arguments!');
           end 
            
        else
           A = varargin{1};
           B = varargin{2};
           if nargin == 3
               dt = varargin{3};
           elseif nargin == 4
               c = varargin{3};
               dt = varargin{4};
           elseif nargin == 5
               c = varargin{3};
               C = varargin{4};
               dt = varargin{5};
           elseif nargin == 6
               c = varargin{3};
               C = varargin{4};
               D = varargin{5};
               dt = varargin{6};
           elseif nargin == 7
               c = varargin{3};
               C = varargin{4};
               D = varargin{5};
               k = varargin{6};
               dt = varargin{7};
           else
               error('Wrong number of input arguments!');
           end 
        end
        
        % get number of states, inputs, and outputs
        states = size(A,1);
        inputs = size(B,2);
        
        outputs = 1;
        if ~isempty(C)
           outputs = size(C,1); 
        end
        
        % instantiate parent class
        obj@contDynamics(name,states,inputs,outputs); 
        
        % assign object properties
        obj.A = A; obj.B = B; obj.c = c; obj.C = C; obj.D = D; obj.k = k;
        obj.dt = dt;
    end
    
    % invoke function observe so that the superclass can access private
    % functions
    function [R, tcomp] = observe(obj,params,options)
        [R, tcomp] = observe@contDynamics(obj,params,options);
    end
end
end

%------------- END OF CODE --------------
