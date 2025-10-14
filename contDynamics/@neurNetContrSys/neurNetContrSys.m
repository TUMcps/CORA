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
%    W1 = [ 0.318 -0.056 ; 0.163 -0.841 ; 1.159 -0.155 ; 0.128 1.189 ; 1.039 -0.415 ; -0.117 -0.659 ; -0.648 0.984 ; -0.038 0.101 ; -0.199 -0.363 ; 0.882 -0.739 ]; 
%    b1 = [ -0.448 ; 0.714 ; -1.316 ; 0.627 ; -1.331 ; -0.131 ; -1.827 ; 0.536 ; 0.693 ; -0.688 ];
%    layers{1} = nnLinearLayer(W1, b1);
%    layers{2} = nnSigmoidLayer();
%    W2 = [ -1.480 -0.317 -1.457 0.487 1.057 0.302 -0.878 0.536 -0.598 -0.690 ; 0.153 0.319 1.322 -1.567 -0.588 0.219 -2.861 -0.342 0.419 0.671 ]; 
%    b2 = [ 0.151 ; -0.991 ];
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

    sys;    % system dynamics
    nn;     % neural network controller
    dt;     % sampling time

end

methods

    % class constructor
    function obj = neurNetContrSys(varargin)

        % 0. check number of input arguments
        assertNarginConstructor([0,3],nargin);

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
        obj@contDynamics(sysCL.name,sysOL.nrOfDims,max(1,sysOL.nrOfInputs-nn.neurons_out),0);
        obj.sys = sysCL; obj.nn = nn; obj.dt = dt;

    end
end
end


% Auxiliary functions -----------------------------------------------------

function [sys,nn,dt] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

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
        if sys.nrOfDims ~= nn.neurons_in
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

    n = sys.nrOfDims; m = nn.neurons_out;
    % instantiate closed-loop system
    if isa(sys, 'nonlinearSys')

        f = @(x,u) [sys.mFile(x(1:n), [x(n+1:n+m); u]); zeros(m, 1)];
        name = [sys.name, 'Controlled'];
        sys = nonlinearSys(name, f, n+m, max(1, sys.nrOfInputs-m));

    elseif isa(sys, 'nonlinParamSys')

        f = @(x,u,p) [sys.mFile(x(1:n), [x(n+1:n+m); u], p); zeros(m, 1)];
        name = [sys.name, 'Controlled'];
        sys = nonlinParamSys(name, f, n+m, max(1, sys.nrOfInputs-m));

    end

end

% ------------------------------ END OF CODE ------------------------------
