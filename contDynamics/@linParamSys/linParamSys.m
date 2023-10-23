classdef linParamSys < contDynamics
% linParamSys class (linear parametric system)
%
% Syntax:
%    obj = linParamSys(A,B)
%    obj = linParamSys(A,B,type)
%    obj = linParamSys(name,A,B)
%    obj = linParamSys(name,A,B,type)
%
% Inputs:
%    name - name of the system
%    A - system matrix
%    B - input matrix
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
% See also: ---

% Authors:       Matthias Althoff
% Written:       23-September-2010
% Last update:   01-November-2017 (constant and varying parameters)
% Last revision: 18-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------
  

properties (SetAccess = private, GetAccess = public)

    A = 1;                          % state matrix
    B = 0;                          % input matrix
    constParam logical = true;      % constant or time-varying parameters
    
    % ...
    stepSize = 1;
    taylorTerms = [];
    mappingMatrixSet = [];
    power = [];
    E = [];
    F = [];
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

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'linParamSys')
%             obj = varargin{1}; return
%         end

        % 2. parse input arguments: varargin -> vars
        [name,A,B,type] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,A,B,type,nargin);

        % 4. compute number of states and inputs
        [name,A,B,type,states,inputs] = aux_computeProperties(name,A,B,type);
        
        % 5. instantiate parent class, assign properties
        obj@contDynamics(name,states,inputs,0); 
        obj.A = A; obj.B = B; obj.constParam = true;
        if strcmp(type,'varParam')
            obj.constParam = false;
        end

    end
end
end


% Auxiliary functions -----------------------------------------------------

function [name,A,B,type] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 4
        throw(CORAerror('CORA:tooManyInputArgs',4));
    end

    % default name and type
    def_name = 'linParamSys';
    type = 'constParam';
    % init properties
    A = []; B = [];

    % no input arguments
    if nargin == 0
        name = def_name;
        return
    end

    % parse depending on whether first input argument is the name
    if ischar(varargin{1})
        % first input argument: name
        [name,A,B,type] = setDefaultValues({def_name,A,B,type},varargin);

    else
        % set default name
        name = def_name;
        [A,B,type] = setDefaultValues({A,B,type},varargin);

    end
    
end

function aux_checkInputArgs(name,A,B,type,n_in)
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
            {type, 'str', {'constParam','varParam'}}
        }]);

        % check if dimensions fit
        if isnumeric(A)
            states = size(A,1);
        elseif isa(A,'intervalMatrix') || isa(A,'matZonotope')
            states = dim(A,1);
        end
        if (isnumeric(B) && ~isscalar(B) && states ~= size(B,1)) ...
                || ( (isa(B,'intervalMatrix') || isa(B,'matZonotope')) ...
                 && states ~= dim(B,1))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['Column dimension of input matrix B must match '...
                'row/column dimension of state matrix A.']));
        end
        
    end

end

function [name,A,B,type,states,inputs] = aux_computeProperties(name,A,B,type)
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

end

% ------------------------------ END OF CODE ------------------------------
