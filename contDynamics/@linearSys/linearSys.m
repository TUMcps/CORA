classdef linearSys < contDynamics
% linearSys - object constructor for linearSys objects
%
% Description:
%    Generates a linear system object according to the following first
%    order differential equations:
%       x'(t) = A x(t) + B u(t) + c + w(t)
%       y(t)  = C x(t) + D u(t) + k + v(t)
%
% Syntax:
%    obj = linearSys()
%    obj = linearSys(A)
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
%    obj - generated linearSys object
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
%               14-December-2022 (TL, property check in inputArgsCheck)
%               15-January-2023 (MW, allow 0 and 1 input arguments)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    A     % system matrix: n x n
    B     % input matrix: n x m
    c     % constant input: n x 1
    C     % output matrix: q x n
    D     % throughput matrix: q x m
    k     % output offset: q x 1
    taylor = [];
    krylov = [];
end

methods
    
    % class constructor
    function obj = linearSys(varargin)
        
        if nargin > 7
            throw(CORAerror('CORA:tooManyInputArgs',7));
        end
        
        % parse name, system matrix, input matrix
        if ~isempty(varargin) && ischar(varargin{1})
            name = varargin{1};
            varargin = varargin(2:end);
        else
            name = 'linearSys'; % default name
        end
        [A,B] = setDefaultValues({[],[]},varargin);
        varargin = varargin(3:end);
        
        % number of states, inputs, and outputs
        states = size(A,1);
        inputs = states;
        outputs = states;
        
        % number of inputs
        if ~isempty(A) && isempty(B)
            B = zeros(states,1);
        end
        if ~isscalar(B)
            inputs = size(B,2);
        end

        % for c, D, and k: overwrite empty entries by default zeros
        % case C = [] is allowed: yields no output computation in code
        [c, C] = setDefaultValues({zeros(states,1), 1}, varargin);
        if isempty(c)
            c = zeros(states, 1);
        end
        varargin = varargin(3:end);
        if ~isempty(C) && ~isscalar(C)
            outputs = size(C,1);
        end

        [D, k] = setDefaultValues( ...
            {zeros(outputs,inputs), zeros(outputs,1)}, varargin);

        inputArgsCheck({ ...
            {A, 'att', 'numeric', 'matrix'}
            {B, 'att', 'numeric', 'matrix'}
            {c, 'att', 'numeric', 'matrix'}
            {C, 'att', 'numeric', 'matrix'}
            {D, 'att', 'numeric', 'matrix'}
            {k, 'att', 'numeric', 'column'}
        });
        
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