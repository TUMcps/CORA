classdef confPredictor
% confPredictor - basic class for set-based predictors
%
% Syntax:
%    sys = confPredictor(fun)
%    sys = confPredictor(fun,type)
%    sys = confPredictor(fun,type,task)
%    sys = confPredictor(fun,type,task,name)
%
% Inputs:
%    fun - nonlinearSysDT object, neural netowrk, or function handle for
%            the point predictor
%    type - predictor type ('zonotope', 'interval', or 'conformalAdd')
%    task - prediction task ('reg" or 'class')
%    name - system name
%
% Outputs:
%    pred - generated confPredictor object
%
% Example:
%    dim_x = 2;
%    dim_u = 2;
%    dim_y = 2;
%    fun = @(x,u) x.^2 + u;
%    sys = nonlinearSysDT(@(x,u) x, 1, dim_x, dim_u, fun, dim_y, false);
%    pred = confPredictor(sys, 'zonotope');
%
% References: 
%    [1] Laura Luetzow, Michael Eichelbeck, Mykel Kochenderfer, and
%        Matthias Althoff. "Zono-Conformal Prediction: Zonotope-Based 
%        Uncertainty Quantification for Regression and Classification 
%        Tasks," arXiv, 2025.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       25-September-2025 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % general properties
    name;               % name of the system
    nrOfInputs;         % input dimension
    nrOfOutputs;        % output dimension
    nrOfUncertainties   % uncertainty dimension
    sys                 % nonlinearSysDT object
    fun                 % underlying prediction function
    type                % predictor type
    task                % prediction task

    % calibration properties
    testSuite           % calibration data as traejctory objects
    U_init              % initial uncertainty set
    U                   % calibrated uncertainty set
    results             % calibration results
end
    
methods

    % class constructor
    function pred = confPredictor(fun,type,task,name)
        %create the zono-conformal predictor from data

        % 0. check number of input arguments
        assertNarginConstructor(1:4,nargin);

        % 1. set default values and check correctness of input arguments
        if nargin < 2
            type = "zonotope";
            task = "reg";
            name = "";
        elseif nargin < 3
            task = "reg";
            name = "";
        elseif nargin < 4
            name = "";
        end
        aux_checkInputArgs(fun,type,task,name);

        % 2. assign properties
        pred.name = name;
        pred.type = type;
        pred.task = task;

        % create nonlinearSysDT object
        if isa(fun,'neuralNetwork')
            pred.fun = fun;
            % transform to contDynamics object
            pred.sys = nonlinearSysDT(fun);
        elseif isa(fun,'function_handle')
            pred.fun = fun;
            pred.sys = nonlinearSysDT(@(x,u) x,1,@(x,u) func(x,u,pred));
        elseif isa(fun, 'nonlinearSysDT')
            pred.sys = fun;
        end
        pred.nrOfInputs = pred.sys.nrOfDims;
        pred.nrOfUncertainties = pred.sys.nrOfInputs;
        pred.nrOfOutputs = pred.sys.nrOfOutputs;
    end
end

methods (Static = true)
    output = evaluateCoverage(varargin)       
    % theoretical coverage guarantees
end

methods (Access = protected)
    [printOrder] = getPrintSystemInfo(S)
end

% getter & setter ---------------------------------------------------------

methods 
    function pred = setType(pred,type)
        pred.type = type;

        % check if predictor type fits to the set representation 
        if ~isempty(pred.U_init) && ...
                ((strcmp(type,"zonotope") && ~isa(pred.U_init,'zonotope')) || ...
                (strcmp(type,"interval") && ~isa(pred.U_init,'interval')))
            fprintf("Set representation of initial uncertainty set " + ...
                "U_init does not fit the new predictor type. " + ...
                "Property U_init will be set to []. \n")
            pred.U_init = [];
        end
    end    

    function pred = setUinit(pred,U_init)
        if (strcmp(pred.type,"zonotope") && isa(U_init,'zonotope')) || ...
                (strcmp(pred.type,"interval") && isa(U_init,'interval'))
            % correct uncertainty set representation
            pred.U_init = U_init;
        else
            throw(CORAerror("CORA:specialError", "Wrong set " + ...
                "representation for initial uncertainty set."))
        end
    end    

    function pred = setU(pred,U)
        if (strcmp(pred.type,"zonotope") && isa(U,'zonotope')) || ...
                (strcmp(pred.type,"interval") && isa(U,'interval')) || ...
                (strcmp(pred.type,"conformalAdd") && strcmp(pred.task,"reg") &&  isa(U,'interval')) || ...
                (strcmp(pred.type,"conformalAdd") && strcmp(pred.task,"class") &&  isnumeric(U))
            % correct uncertainty set representation
            pred.U = U;
        else
            throw(CORAerror("CORA:specialError", "Wrong set " + ...
                "representation for uncertainty set."))
        end
    end    
end
end


% Auxiliary functions -----------------------------------------------------

function aux_checkInputArgs(fun,type,task,name)

% check fun 
if ~isa(fun,'nonlinearSysDT') && ~isa(fun,'neuralNetwork') && ...
        ~(isa(fun,'function_handle') && nargin(fun) == 2)
    throw(CORAerror('CORA:wrongInputInConstructor',...
        ['Prediction function has to be a nonlinearSysDT object, a ' ...
        'neural network or a function handle with two input arguments.']));
end

% check type 
if ~strcmp(type,'zonotope') && ~strcmp(type,'interval') && ~strcmp(type,'conformalAdd')
    throw(CORAerror('CORA:wrongInputInConstructor',...
        ['Prediction type must be interval or zonotope.']));
end

% check task 
if ~strcmp(task,'reg') && ~strcmp(task,'class')
    throw(CORAerror('CORA:wrongInputInConstructor',...
        ['Prediction task must be reg for regression tasks or ' ...
        'class for classification tasks.']));
end

% check name
if ~ischar(name) && ~isstring(name)
    throw(CORAerror('CORA:wrongInputInConstructor',...
        'System name has to be a char or string array.'));
end    
end

% ------------------------------ END OF CODE ------------------------------
