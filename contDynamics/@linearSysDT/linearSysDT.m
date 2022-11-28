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
%       x(k+1) = A x(k) + B u(k) + c + w(k)
%       y(k) = C x(k) + D u(k) + k + v(k)
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
%    obj - generated linearSysDT object
%
% Example:
%    A = [-0.4 0.6; 0.6 -0.4];
%    B = [0; 1];
%    C = [1 0];
%    dt = 0.4;
%    sys = linearSysDT(A,B,0,C,dt)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      20-Mar-2020 
% Last update:  14-Jun-2021 (MA, invoke observe from superclass)
%               19-November-2021 (MW, default values for C, D, k)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    A {mustBeNumeric} = [];                         % system matrix: n x n
    B {mustBeNumeric} = 1;                          % input matrix: n x m
    c {mustBeNumeric} = [];                     	% constant input: n x 1
    C {mustBeNumeric} = [];                         % output matrix: q x n
    D {mustBeNumeric} = [];                         % throughput matrix: q x m
    k {mustBeNumeric} = [];                         % output offset: q x 1
    dt {mustBeNumeric} = [];                        % sampling time
end

methods
    
    % class constructor
    function obj = linearSysDT(varargin)
        
        % not enough or too many input arguments
        if nargin < 3 || nargin > 8
            throw(CORAerror('CORA:notEnoughInputArgs',3));
        elseif nargin > 8
            throw(CORAerror('CORA:tooManyInputArgs',8));
        end

        % parse name, system matrix, input matrix
        if ischar(varargin{1})
            name = varargin{1};
            A = varargin{2};
            B = varargin{3};
            cnt = 3;
        else
            name = 'linearSysDT'; % default name
            A = varargin{1};
            B = varargin{2};
            cnt = 2;
        end
        % note that cnt is +1 compared to linearSys due to dt-arguments
        
        % get number of states, inputs, and outputs
        states = size(A,1);
        inputs = states;
        outputs = states;
        
        % number of inputs
        if ~isscalar(B)
            inputs = size(B,2);
        end
        
        % for c, D, and k: overwrite empty entries by default zeros
        % case C = [] is allowed: yields no output computation in code
        
        % constant input
        c = zeros(states,1); % default value
        if nargin > cnt + 1
        	cnt = cnt + 1;
            if ~isempty(varargin{cnt})
                c = varargin{cnt};
            end 
        end
        
        % output matrix
        C = 1; % default value
        if nargin > cnt + 1
        	cnt = cnt + 1;
            C = varargin{cnt};
        end
        % compute number of outputs
        if ~isempty(C) && ~isscalar(C)
            outputs = size(C,1); 
        end
        
        % throughput matrix
        D = 0; % default value
        if inputs ~= outputs
            D = zeros(outputs,inputs); % default value
        end
        if nargin > cnt + 1
        	cnt = cnt + 1;
            if ~isempty(varargin{cnt})
                D = varargin{cnt};
            end 
        end
        
        % output offset
        k = zeros(outputs,1); % default value
        if nargin > cnt + 1
        	cnt = cnt + 1;
            if ~isempty(varargin{cnt})
                k = varargin{cnt};
            end
        end
        
        % sampling time (always last input argument)
        dt = varargin{end};
        
        % instantiate parent class
        obj@contDynamics(name,states,inputs,outputs); 
        
        % assign object properties
        obj.A = A; obj.B = B; obj.c = c;
        obj.C = C; obj.D = D; obj.k = k;
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
