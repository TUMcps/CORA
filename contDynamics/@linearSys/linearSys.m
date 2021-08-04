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
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    A = [];             % system matrix
    B = 1;              % input matrix
    c = [];             % constant input
    C = 1;              % output matrix
    D = [];             % throughput matrix
    k = [];             % output offset
    taylor = [];
    krylov = [];
end

methods
    
    % class constructor
    function obj = linearSys(varargin)
        
        c = []; C = 1; D = []; k = [];
        
        % parse input arguments
        if ischar(varargin{1})
            name = varargin{1};
            A = varargin{2};
            B = varargin{3};
            cnt = 3;
        else
            name = 'linearSys';
            A = varargin{1};
            B = varargin{2};
            cnt = 2;
        end
        
        if nargin > cnt
           cnt = cnt + 1;
           c = varargin{cnt}; 
        end
        if nargin > cnt
           cnt = cnt + 1;
           C = varargin{cnt}; 
        end
        if nargin > cnt
           cnt = cnt + 1;
           D = varargin{cnt}; 
        end
        if nargin > cnt
           cnt = cnt + 1;
           k = varargin{cnt}; 
        end
        
        % number of states, inputs, and outputs
        states = size(A,1);
        if ~isscalar(B)
            inputs = size(B,2);
        else
            inputs = states; 
        end
        outputs = 1;
        if ~isempty(C)
           outputs = size(C,1); 
        end
        
        % instantiate parent class
        obj@contDynamics(name,states,inputs,outputs); 
        
        % assign object properties
        obj.A = A; obj.B = B; obj.c = c; obj.C = C; obj.D = D; obj.k = k;
    end
end
end

%------------- END OF CODE --------------