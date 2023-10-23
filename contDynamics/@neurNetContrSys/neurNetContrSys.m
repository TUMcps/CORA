classdef neurNetContrSys < contDynamics
% neurNetContrSys - class that stores neural network controlled systems
%
% Syntax:
%    obj = neurNetContrSys(sysOL,nn,dt)
%
% Inputs:
%    sysOL - dynamics of the uncontrolled system (class: contDynamics)
%    nn - neural network controller (class: neuralNetwork)
%    dt - sampling time
%
% Outputs:
%    obj - generated neurNetContrSys object
%
% Example:
%    % dynamic system
%    f = @(x,u) [x(2) + u(2); (1-x(1)^2)*x(2) - x(1) + u(1)];
%    sysOL = nonlinearSys(f);
%
%    % neural network controller
%    layers = cell(4, 1);
%    W1 = rand(10,2); b1 = rand(10,1);
%    layers{1} = nnLinearLayer(W1, b1);
%    layers{2} = nnSigmoidLayer();
%    W2 = rand(2,10); b2 = rand(2,1);
%    layers{3} = nnLinearLayer(W2, b2);
%    layers{4} = nnSigmoidLayer();
%    nn = neuralNetwork(layers);
%
%    % neural network controlled system
%    dt = 0.01;
%    sys = neurNetContrSys(sysOL,nn,dt);
%
% Reference:
%   [1] Kochdumper, Niklas, et al. "Open-and Closed-Loop Neural Network
%       Verification using Polynomial Zonotopes." arXiv preprint 
%       arXiv:2207.02715 (2022).
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neurNetContrSys

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   23-November-2022 (TL, polish)
%                14-December-2022 (TL, property check in inputArgsCheck)
% Last revision: 18-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    sys; % system dynamics
    nn; % neural network controller
    dt; % sampling time
end

methods

    % class constructor
    function obj = neurNetContrSys(varargin)

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'neurNetContrSys')
%             obj = varargin{1}; return
%         end

        % 2. parse input arguments: varargin -> vars
        [sysOL,nn,dt] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(sysOL,nn,dt,nargin);

        % 4. instantiate closed-loop system, convert old neuralNetwork
        [sysCL,nn,dt] = aux_computeProperties(sysOL,nn,dt);
        
        % 5. instantiate parent class, assign properties
        obj@contDynamics(sysCL.name,sysOL.dim,max(1,sysOL.nrOfInputs-nn.neurons_out),0);
        obj.sys = sysCL; obj.nn = nn; obj.dt = dt;

    end
end
end


% Auxiliary functions -----------------------------------------------------

function [sys,nn,dt] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin ~= 0 && nargin < 3
        throw(CORAerror('CORA:notEnoughInputArgs',3));
    elseif nargin > 3
        throw(CORAerror('CORA:tooManyInputArgs',3));
    end

    % default values
    sys = contDynamics(); nn = []; dt = 0;

    % no input arguments
    if nargin == 0
        return
    end

    % parse user-provided input arguments
    [sys,nn,dt] = setDefaultValues({sys,nn,dt},varargin);
    
end

function aux_checkInputArgs(sys,nn,dt,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % check data types
        inputArgsCheck({ ...
            {sys,'att',{'nonlinearSys','nonlinParamSys'}},...
            {nn,'att',{'neuralNetwork'}},...
            {dt,'att','numeric',{'scalar','positive'}},...
        })

        % check if dimensions fit
        if sys.dim ~= nn.neurons_in
            throw(CORAerror('CORA:wrongInputInConstructor', ...
               'Dimension of sys and input of nn should match.'));
        end
        if sys.nrOfInputs < nn.neurons_out
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['Dimensions of open-loop system and neural network', ...
                'are not consistent!']));
        end
        
    end

end

function [sys,nn,dt] = aux_computeProperties(sys,nn,dt)
% compute properties of neurNetContrSys object

    n = sys.dim; m = nn.neurons_out;
    % instantiate closed-loop system
    if isa(sys, 'nonlinearSys')

        f = @(x, u) [sys.mFile(x(1:n), [x(n+1:n+m); u]); zeros(m, 1)];
        name = [sys.name, 'Controlled'];
        sys = nonlinearSys(name, f, n+m, max(1, sys.nrOfInputs-m));

    elseif isa(sys, 'nonlinParamSys')

        f = @(x, u, p) [sys.mFile(x(1:n), [x(n+1:n+m); u], p); zeros(m, 1)];
        name = [sys.name, 'Controlled'];
        sys = nonlinParamSys(name, f, n+m, max(1, sys.nrOfInputs-m));

    end

end

% ------------------------------ END OF CODE ------------------------------
