classdef linProbSys < contDynamics
% linProbSys - class constructor for linear probabilistic systems
%
% Syntax:
%    obj = linProbSys(A,B)
%    obj = linProbSys(A,B,C)
%    obj = linProbSys(name,A,B)
%    obj = linProbSys(name,A,B,C)
%
% Inputs:
%    name - name of system
%    A - state matrix
%    B - input matrix
%    C - noise matrix
%
% Outputs:
%    obj - generated linProbSys object
%
% Example:
%    A = [-1 -4; 4 -1];
%    B = eye(2);
%    C = 0.7*eye(2);
%    sys = linProbSys(A,B,C);
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
%               19-June-2022 (MW, update syntax)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    A {mustBeNumeric} = []; % state matrix
    B {mustBeNumeric} = []; % input matrix
    C {mustBeNumeric} = []; % noise matrix
    taylor = [];
end

methods
    % class constructor
    function obj = linProbSys(varargin)

        % check number of input arguments
        if nargin < 2
            throw(CORAerror('CORA:notEnoughInputArgs',2));
        elseif nargin > 4
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end

        % parse name, set index for other input arguments
        if ischar(varargin{1})
            name = varargin{1};
            cnt = 2;
        else
            % default name
            name = 'linProbSys';
            cnt = 1;
        end

        % parse state matrix
        A = varargin{cnt};
        if nargin == 2 && ischar(varargin{1})
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['If the first input argument is the name of the system, '...
                'there have to be either 3 or 4 input arguments.']));
        end

        % parse input matrix
        B = varargin{cnt+1};

        % parse noise matrix
        if nargin == cnt+2
            C = varargin{cnt+2};
        end

        % instantiate parent class (never any outputs)
        obj@contDynamics(name,size(A,1),size(B,2),0);

        % assign properties
        obj.A = A;
        obj.B = B;
        if nargin == cnt+2
            obj.C = C;
        end

    end
end
end

%------------- END OF CODE --------------