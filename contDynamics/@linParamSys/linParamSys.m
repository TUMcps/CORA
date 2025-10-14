classdef linParamSys < contDynamics
% linParamSys class (linear parametric system)
%
% Syntax:
%    obj = linParamSys(A)
%    obj = linParamSys(A,B)
%    obj = linParamSys(A,B,c)
%    obj = linParamSys(A,type)
%    obj = linParamSys(A,B,type)
%    obj = linParamSys(A,B,c,type)
%    obj = linParamSys(name,A,type)
%    obj = linParamSys(name,A,B,type)
%    obj = linParamSys(name,A,B,c,type)
%
% Inputs:
%    name - name of the system
%    A - system matrix
%    B - input matrix
%    c - constant input
%    type - constant/time-varying parameters
%           - 'constParam' (constant parameters, default)
%           - 'varParam' (time-varying parameters)
%
% Outputs:
%    obj - generated linParamSys object
%
% Example:
%    Ac = [-2 0; 1.5 -3]; Aw = [0 0; 0.5 0];
%    A = intervalMatrix(Ac,Aw);
%    B = [1; 1];
%    sys = linParamSys(A,B,'varParam')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       23-September-2010
% Last update:   01-November-2017 (constant and varying parameters)
%                10-September-2025 (MP, Allow autonomous systems)
%                10-September-2025 (MP, Add constant input c)
% Last revision: 18-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------
  

properties (SetAccess = private, GetAccess = public)

    A = 1;                          % state matrix
    B = 0;                          % input matrix
    c = 0;                          % constant input
    type = '';
    constParam logical = true;      % constant or time-varying parameters
    
    % TODO: check if these are still necessary...
    stepSize = 1;
    taylorTerms = [];
    mappingMatrixSet = [];
    power = [];
    E = []; % note: not disturbance matrix
    F = []; % note: not noise matrix
    inputF = [];
    inputCorr = [];
    Rinput = [];
    Rtrans = [];
    RV = [];
    sampleMatrix = [];
end
    
methods
    
    % class constructor
    function obj = linParamSys(varargin)

        % 0. check number of input arguments
        assertNarginConstructor(0:5,nargin);

        % 1. parse input arguments: varargin -> vars
        [name,A,B,c,type] = aux_parseInputArgs(varargin{:});

        % 2. check correctness of input arguments
        aux_checkInputArgs(name,A,B,c,type,nargin);

        % 3. compute number of states and inputs
        [name,A,B,c,type,states,inputs] = aux_computeProperties(name,A,B,c,type);
        
        % 4. instantiate parent class, assign properties
        obj@contDynamics(name,states,inputs,0); 
        obj.A = A; obj.B = B; obj.c = c; obj.type = type; obj.constParam = true;
        if strcmp(type,'varParam')
            obj.constParam = false;
        end

    end
end

methods (Access = protected)
    [printOrder] = getPrintSystemInfo(S)
end

end


% Auxiliary functions -----------------------------------------------------

function [name,A,B,c,type] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % default name and type
    def_name = 'linParamSys'; type = 'constParam'; A = []; B = []; c = [];

    def_type = 'constParam';
    % no input arguments
    if nargin == 0
        name = def_name;
        return
    end

    
    % parse depending on whether first input argument is the name
    if ischar(varargin{1})
        % edge case: if size(varargin) == 3/4, it is ambigious wether we have
        % [[name,]A,B,c] or [[name,]A,B,type]
        if size(varargin,2) == 4 && ischar(varargin{4})
            varargin{5} = varargin{4};
            varargin{4} = c;
        end

        % first input argument: name
        [name,A,B,c,type] = setDefaultValues({def_name,A,B,c,type},varargin);

    else
        % same as above
        if size(varargin,2) == 3 && ischar(varargin{3})
            varargin{4} = varargin{3};
            varargin{3} = c;
        end
        % set default name
        name = def_name;
        [A,B,c,type] = setDefaultValues({A,B,c,type},varargin);

    end
    
end

function aux_checkInputArgs(name,A,B,c,type,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % ensure that values have correct data type
        inputChecks = {};
        if strcmp(name,'linParamSys')
            % default name (unless explicitly chosen by user, we have A as
            % first input argument)
            inputChecks = {{name,'att',{'char','string'}}};
        end
        inputArgsCheck([inputChecks;
            {{A, 'att', {'numeric','intervalMatrix','matZonotope'}}
            {B, 'att', {'numeric','intervalMatrix','matZonotope'}}
            {c, 'att', {'numeric','interval','zonotope'}}
            {type, 'str', {'constParam','varParam'}}
        }]);

        % check if dimensions fit
        if isnumeric(A)
            states = size(A,1);
        elseif isa(A,'intervalMatrix') || isa(A,'matZonotope')
            states = dim(A,1);
        end
        if ~isempty(B) && ((isnumeric(B) && ~isscalar(B) && states ~= size(B,1)) ...
                || ( (isa(B,'intervalMatrix') || isa(B,'matZonotope')) ...
                 && states ~= dim(B,1)))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['Column dimension of input matrix B must match '...
                'row/column dimension of state matrix A.']));
        end 
        if ~isempty(c) && ~isvector(c)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Offset c must be a vector.'));
        end
        % check dimension of offset for state
        if ~isempty(c) && states ~= length(c)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Length of offset c must match row/column dimension of state matrix A.'));
        end

    end

end

function [name,A,B,c,type,states,inputs] = aux_computeProperties(name,A,B,c,type)
% assign zero vectors/matrices for [] values (from user or by default)
% compute number of states and inputs
% difficulty: MATLAB can use 1 instead of eye(n) for multiplication

    % number of states
    if isnumeric(A)
        states = size(A,1);
    elseif isa(A,'intervalMatrix') || isa(A,'matZonotope')
        states = dim(A,1);
    end

    % input matrix
    if ~isempty(A) && isempty(B)
        B = zeros(states,1);
    end

    % number of inputs
    if isscalar(B)
        inputs = states;
    elseif isa(B,'matZonotope') || isa(B,'intervalMatrix')
        inputs = dim(B,2);
    else
        inputs = size(B,2);
    end

    % constant offset
    if isempty(c)
        c = zeros(states,1);
    end
end

% ------------------------------ END OF CODE ------------------------------
