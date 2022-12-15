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

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      17-September-2021
%               23-November-2022 (TL, polish)
%               14-December-2022 (TL, property check in inputArgsCheck)
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)

    sys; % system dynamics
    nn; % neural network controller
    dt; % sampling time
end

methods

    % class constructor
    function obj = neurNetContrSys(sys,nn,dt)

        % check user input
        if nargin < 3
            throw(CORAerror('CORA:notEnoughInputArgs',3));
        elseif nargin > 3
            throw(CORAerror('CORA:tooManyInputArgs',3));
        end
        if isa(nn, "neuralNetworkOld")
            nn = neuralNetwork.getFromOldNeuralNetwork(nn);
        end

        inputArgsCheck({ ...
            {sys, 'att', 'contDynamics'}; ...
            {nn, 'att', 'neuralNetwork'}; ...
            {dt, 'att', 'numeric', {'finite', 'nonnegative', 'scalar'}};
        });

        % check dimensions
        n = sys.dim;
        m = nn.neurons_out;

        if n ~= nn.neurons_in
            throw(CORAerror('CORA:wrongInputInConstructor', ...
               'Dimension of sys and input of nn should match.'));
        end
        if sys.nrOfInputs < m
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['Dimensions of open-loop system and neural network', ...
                'are not consistent!']));
        end

        % construct closed-loop system
        if isa(sys, 'nonlinearSys')

            f = @(x, u) [sys.mFile(x(1:n), [x(n+1:n+m); u]); zeros(m, 1)];
            name = [sys.name, 'Controlled'];
            sys = nonlinearSys(name, f, n+m, max(1, sys.nrOfInputs-m));

        elseif isa(sys, 'nonlinParamSys')

            f = @(x, u, p) [sys.mFile(x(1:n), [x(n+1:n+m); u], p); zeros(m, 1)];
            name = [sys.name, 'Controlled'];
            sys = nonlinParamSys(name, f, n+m, max(1, sys.nrOfInputs-m));

        else
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Only nonlinearSysa and nonlinParamSys supported.'));
        end

        % generate parent object
        obj@contDynamics(sys.name, n, max(1, sys.nrOfInputs-m), 1);

        % assign object properties
        obj.sys = sys;
        obj.nn = nn;
        obj.dt = dt;
    end
end
end

%------------- END OF CODE --------------