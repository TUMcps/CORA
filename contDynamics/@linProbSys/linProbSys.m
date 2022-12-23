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
%               14-December-2022 (TL, property check in inputArgsCheck)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    A; % state matrix
    B; % input matrix
    C; % noise matrix
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
            varargin = varargin(2:end);
        else
            % default name
            name = 'linProbSys';
        end

        if length(varargin) < 3
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['If the first input argument is the name of the system, '...
                'there have to be either 3 or 4 input arguments.']));
        end

        [A, B, C] = setDefaultValues({[], [], []}, varargin);

        inputArgsCheck({ ...
            {A, 'att', 'numeric'}
            {B, 'att', 'numeric'}
            {C, 'att', 'numeric'}
        });

        % instantiate parent class (never any outputs)
        obj@contDynamics(name,size(A,1),size(B,2),0);

        % assign properties
        obj.A = A;
        obj.B = B;
        obj.C = C;
    end
end
end

%------------- END OF CODE --------------