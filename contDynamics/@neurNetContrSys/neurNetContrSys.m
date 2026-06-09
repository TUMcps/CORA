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
%                20-February-2026 (TL, add linearSys/linearSysDT/nonlinearSysDT/linParamSys)
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
            {sys,'att',{'linearSys','linearSysDT','linParamSys','nonlinearSys','nonlinearSysDT','nonlinParamSys'}},...
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
    name = [sys.name, 'Controlled'];
    n_ext = max(1, sys.nrOfInputs - m);

    % instantiate closed-loop system
    if isa(sys, 'nonlinearSys') % nonlinear ---

        f = @(x,u) [sys.mFile(x(1:n), [x(n+1:n+m); u]); zeros(m, 1)];
        sys = nonlinearSys(name, f, n+m, n_ext);

    elseif isa(sys, 'nonlinearSysDT')

        f = @(x,u) [sys.mFile(x(1:n), [x(n+1:n+m); u]); zeros(m, 1)];
        sys = nonlinearSysDT(name, f, dt, n+m, n_ext);

    elseif isa(sys, 'nonlinParamSys')

        f = @(x,u,p) [sys.mFile(x(1:n), [x(n+1:n+m); u], p); zeros(m, 1)];
        sys = nonlinParamSys(name, f, n+m, n_ext);

    elseif isa(sys, 'linearSys') % linear ---

        [A_aug, B_aug, c_aug] = aux_augmentLinear(sys, n, m, n_ext);
        sys = linearSys(name, A_aug, B_aug, c_aug);

    elseif isa(sys, 'linearSysDT')

        [A_aug, B_aug, c_aug] = aux_augmentLinear(sys, n, m, n_ext);
        sys = linearSysDT(name, A_aug, B_aug, c_aug, dt);

    elseif isa(sys, 'linParamSys')

        [A_aug, B_aug, c_aug] = aux_augmentLinParam(sys, n, m, n_ext);
        sys = linParamSys(name, A_aug, B_aug, c_aug, sys.type);

    end

end

function [A_aug, B_aug, c_aug] = aux_augmentLinear(sys, n, m, n_ext)
% build augmented matrices for linear closed-loop system
%   augmented state: [x; u_nn], where u_nn = x(n+1:n+m) are the NN outputs
%   x'     = A*x + B_nn*u_nn + B_ext*u_ext + c  =>  [A, B_nn; 0, 0] * x_aug + [B_ext; 0] * u_ext
%   u_nn'  = 0  (reset by NN at each sampling step)

    % expand scalar B to full matrix (scalar B means identity-like effect)
    B = sys.B;
    if isscalar(B)
        B = B * eye(n);
    end

    A_aug = [sys.A, B(:, 1:m); zeros(m, n+m)];

    if sys.nrOfInputs > m
        B_aug = [B(:, m+1:end); zeros(m, n_ext)];
    else
        B_aug = zeros(n+m, 1);
    end

    c_aug = [sys.c; zeros(m, 1)];

end

function [A_aug, B_aug, c_aug] = aux_augmentLinParam(sys, n, m, n_ext)
% build augmented matrices for linParamSys closed-loop system
%   same structure as aux_augmentLinear; A may be intervalMatrix/matZonotope,
%   B is assumed numeric (uncertainty typically lives in A only)

    B = sys.B;
    if isscalar(B)
        B = B * eye(n);
    end

    B_nn = B(:, 1:m);

    if isa(sys.A, 'intervalMatrix')
        % intervalMatrix has no horzcat/vertcat; work through numeric bounds
        A_inf = infimum(sys.A.int);
        A_sup = supremum(sys.A.int);
        A_aug_inf = [A_inf, B_nn; zeros(m, n+m)];
        A_aug_sup = [A_sup, B_nn; zeros(m, n+m)];
        A_aug = intervalMatrix((A_aug_inf+A_aug_sup)/2, (A_aug_sup-A_aug_inf)/2);
    elseif isa(sys.A, 'matZonotope')
        % matZonotope: augment center and each generator slice (n x n x h)
        C_aug = [sys.A.C, B_nn; zeros(m, n+m)];
        h = size(sys.A.G, 3);
        G_aug = zeros(n+m, n+m, h);
        for i = 1:h
            G_aug(:,:,i) = [sys.A.G(:,:,i), zeros(n,m); zeros(m,n+m)];
        end
        A_aug = matZonotope(C_aug, G_aug);
    else
        A_aug = [sys.A, B_nn; zeros(m, n+m)];
    end

    if sys.nrOfInputs > m
        B_aug = [B(:, m+1:end); zeros(m, n_ext)];
    else
        B_aug = zeros(n+m, 1);
    end

    c_aug = [sys.c; zeros(m, 1)];

end

% ------------------------------ END OF CODE ------------------------------
