classdef linearSys < contDynamics
% linearSys class (linSys: linear system)
%
% Syntax:  
%    obj = linearSys(A,B)
%    obj = linearSys(A,B,c)
%    obj = linearSys(A,B,c,C)
%    obj = linearSys(A,B,c,C,D)
%    obj = linearSys(A,B,c,C,D,k)
%    obj = linearSys(name,A,B)
%    obj = linearSys(name,A,B,c)
%    obj = linearSys(name,A,B,c,C)
%    obj = linearSys(name,A,B,c,C,D)
%    obj = linearSys(name,A,B,c,C,D,k)
%
% Description:
%    Generates a linear system object according to the following first
%    order differential equations:
%       x'(t) = A x(t) + B u(t) + c + w(t)
%       y(t)  = C x(t) + D u(t) + k + v(t)
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
%    obj - Generated Object
%
% Example:
%    A = [-2 0; 1 -3];
%    B = [1; 1];
%    C = [1 0];
%
%    sys = linearSys(A,B,[],C)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSysDT

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      23-January-2007 
% Last update:  30-April-2007
%               04-August-2016 (changed to new OO format)
%               01-November-2017 (constant input added)
%               20-March-2018 (NK, output equation parameter added)
%               07-November-2018 (MA, default values for B and C changed)
%               04-March-2019 (MA, default IDs for inputs and outputs)
%               20-May-2020 (NK, name now optional)
%               19-November-2021 (MW, default values for c, D, k)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    A {mustBeNumeric} = [];     % system matrix: n x n
    B {mustBeNumeric} = 1;      % input matrix: n x m
    c {mustBeNumeric} = [];     % constant input: n x 1
    C {mustBeNumeric} = [];     % output matrix: q x n
    D {mustBeNumeric} = [];     % throughput matrix: q x m
    k {mustBeNumeric} = [];     % output offset: q x 1
    taylor = [];
    krylov = [];
end

methods
    
    % class constructor
    function obj = linearSys(varargin)
        
        if nargin < 2
            % not enough input arguments
            throw(CORAerror('CORA:notEnoughInputArgs',2));
        elseif nargin > 7
            % too many input arguments
            throw(CORAerror('CORA:tooManyInputArgs',7));
        end
        
        % parse name, system matrix, input matrix
        if ischar(varargin{1})
            name = varargin{1};
            A = varargin{2};
            B = varargin{3};
            cnt = 3;
        else
            name = 'linearSys'; % default name
            A = varargin{1};
            B = varargin{2};
            cnt = 2;
        end
        
        % number of states, inputs, and outputs
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
        if nargin > cnt
        	cnt = cnt + 1;
            if ~isempty(varargin{cnt})
                c = varargin{cnt};
            end
        end
        
        % output matrix
        C = 1; % default value
        if nargin > cnt
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
        if nargin > cnt
            cnt = cnt + 1;
            if ~isempty(varargin{cnt})
                D = varargin{cnt};
            end
        end
        
        % output offset
        k = zeros(outputs,1); % default value
        if nargin > cnt
        	cnt = cnt + 1;
            if ~isempty(varargin{cnt})
                k = varargin{cnt};
            end
        end
        
        % instantiate parent class
        obj@contDynamics(name,states,inputs,outputs); 
        
        % assign object properties
        obj.A = A; obj.B = B; obj.c = c;
        obj.C = C; obj.D = D; obj.k = k;
    end
end

methods (Static = true)
    linSys = generateRandom(varargin) % generates random linear system
end

end

%------------- END OF CODE --------------