classdef linearARX < contDynamics
% linearARX - object constructor for linear discrete-time ARX systems
%
% Syntax:
%    linARX = linearARX(A_bar,B_bar,dt)
%    linARX = linearARX(name,A_bar,B_bar,dt)
%
% Description:
%    Generates a discrete-time linear ARX object according to the 
%    following equation:
%       y(k) =  \sum_{i=1}^p A_bar{i} y(k-i) + 
%               \sum_{i=1}^{p+1} B_bar{i} u(k-i+1)
%
% Inputs:
%    name   - name of system
%    A_bar  - output parameters
%    B_bar  - input parameters
%    dt     - sampling time
%
% Outputs:
%    linARX - linearARX object
%
% Example:
%    dt = 0.1;
%    A_bar = {[-0.4 0.6; 0.6 -0.4];[0.1 0; 0.2 -0.5]};
%    B_bar = {[0; 0];[0.3; -0.7];[0.1; 0]};
%    sys = linearARX(A_bar,B_bar,dt)
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       02-February-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    A_bar       % output parameters
    B_bar       % input parameters
    dt          % sampling time
    tvp         % use of normal parametrization or time-varying parameters
    conv_tvp    % struct containing the parameters for conversion to tvp
    A_tilde     % time-varying output parameters
    B_tilde     % time-varying input parameters
    n_p         % number of past time steps
end

methods

    % class constructor
    function linARX = linearARX(varargin)

        % 0. check number of input arguments
        assertNarginConstructor(3:4,nargin);

        % 1. copy constructor: not allowed due to obj@contDynamics below
%         if nargin == 1 && isa(varargin{1},'linearARX')
%             obj = varargin{1}; return
%         end

        % 2. parse input arguments: varargin -> vars
        [name,A_bar,B_bar,dt] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,A_bar,B_bar,dt,nargin);

        % 4. compute number of states, inputs, and outputs
        [name,A_bar,B_bar,dt,n_p,inputs,outputs] = ...
            aux_computeProperties(name,A_bar,B_bar,dt);
        
        % 5a. instantiate parent class
        linARX@contDynamics(name,0,inputs,outputs);

        % 5b. assign object properties
        linARX.A_bar = A_bar;
        linARX.B_bar = B_bar;
        linARX.dt = dt;
        linARX.tvp = false;
        linARX.n_p = n_p;
    end

    % invoke function observe so that the superclass can access private
    % functions
    function [R, tcomp] = observe(linARX,params,options)
        [R, tcomp] = observe@contDynamics(linARX,params,options);
    end

    function linARX = setTVP(linARX)
        % set tvp attribute to true and compute the parameters required for 
        % the computation of the time-varying parameters A_tilde and B_tilde

        linARX.tvp = true;
        if ~isfield(linARX, "conv_tvp")
            p = linARX.n_p;
            dim_y = linARX.nrOfOutputs;

            % build augmented matrices
            linARX.conv_tvp.A_ext = [zeros(dim_y*(p-1), dim_y) eye(dim_y * (p-1)); zeros(dim_y, dim_y*p)];
            linARX.conv_tvp.B_ext = cell(p+1,1);
            linARX.conv_tvp.E = [zeros(dim_y, dim_y * (p-1)) eye(dim_y)];
            for i = 0:p
                if i < p
                    col_i = i*dim_y+1:(i+1)*dim_y;
                    linARX.conv_tvp.A_ext(end-dim_y+1:end,col_i) = linARX.A_bar{p-i};
                end
                linARX.conv_tvp.B_ext{i+1} = linARX.conv_tvp.E'*linARX.B_bar{i+1};
            end
        end
    end

    function linARX = negateTVP(linARX)
        % set tvp attribute to false
        linARX.tvp = false;
    end

    function linARX = computeTVP(linARX, k, iVec)
        % compute time-varying parameters A_tilde and B_tilde for the time 
        % point k and the indexes in iVec

        p = linARX.n_p;
        k_plus = k+p-1;
        dim_u = linARX.nrOfInputs;
        dim_y = linARX.nrOfOutputs;
        if ~linARX.tvp
            linARX.setTVP;
        end

        % check if A_tilde
        if isfield(linARX, "A_tilde")
            k_last = length(linARX.A_tilde,1);
        else
            k_last = 0;
            linARX.A_tilde = cell(k+1,1);
            linARX.B_tilde = cell(max(iVec)+1,k+1);
        end

        % check if valid construction
        if max(iVec) > k_plus || min(iVec) < 0
            throw(CORAerror("CORA:notSupported",...
                "Invalid. Index i must be in 0<=i<=k_plus."));
        end
        if k < p
            throw(CORAerror("CORA:notSupported",...
                "Invalid. Time step k must be in k>=p."));
        end
        
        if k == 1
            linARX.A_tilde{k+1} = linARX.conv_tvp.A_ext;
            for i = iVec
                linARX.B_tilde{i+1,k+1} = linARX.conv_tvp.B_ext{i+1};
            end
        else 
            if k_last == 0 
                linARX.A_tilde{k+1} = linARX.conv_tvp.A_ext^(k); % [1, Eq. (9a)] 
            else
                linARX.A_tilde{k+1} = linARX.A_tilde{k_last} * ...
                    linARX.conv_tvp.A_ext^(k-k_last); % [1, Eq. (14a)] 
            end

            % compute parameters from scratch (see [1, Eq. (9b)])
            for i = iVec
                j_min = max(0,i-p);
                j_max = min(k-1, i);
                linARX.B_tilde{i+1,k+1} = zeros(p*dim_y, dim_u);
                Aj = linARX.conv_tvp.A_ext^j_min;
                for j=j_min:j_max
                    if j > j_min 
                        Aj = Aj * linARX.conv_tvp.A_ext;
                    end
                    linARX.B_tilde{i+1,k+1} = linARX.B_tilde{i+1,k+1} + Aj * linARX.conv_tvp.B_ext{i-j+1};
                end
            end
        end
    end

    % update system dynamics for the new augmented input [u; w] where w is
    % the process noise acting on all states 
    function linARX = augment_u_with_w(linARX)
        throw(CORAerror('CORA:notSupported','Not implemented for linearARX.'))
    end

    % update system dynamics for the new augmented input [u; v] where v is
    % the measurement noise acting on all outputs 
    function linARX = augment_u_with_v(linARX)
        dim_y = size(linARX.A_bar,1);
        linARX.B_bar{1} = [linARX.B_bar{1} eye(dim_y)];
        for j=2:length(linARX.B_bar)
            linARX.B_bar{j} = [linARX.B_bar{j} zeros(dim_y)];
        end
        linARX.nrOfInputs = linARX.nrOfInputs + dim_y;
    end
end

methods (Access = protected)
    [printOrder] = getPrintSystemInfo(S)
end

end


% Auxiliary functions -----------------------------------------------------

function [name,A_bar,B_bar,dt] = aux_parseInputArgs(varargin)
% parse name, output parameters, input parameters

    def_name = 'linearARX'; A_bar = []; B_bar = []; dt = [];

    % parse depending on whether first input argument is the name
    if ischar(varargin{1}) || isa(varargin{1},'string')
        [name,A_bar,B_bar,dt] = setDefaultValues({def_name,A_bar,B_bar,dt},varargin);
    else
        name = def_name;
        [A_bar,B_bar,dt] = setDefaultValues({A_bar,B_bar,dt},varargin);
    end

end

function aux_checkInputArgs(name,A_bar,B_bar,dt,n_in)
% check correctness of input arguments

% only check if macro set to true
if CHECKS_ENABLED && n_in > 0

    if strcmp(name,'linearARX')
        inputArgsCheck({ ...
                {name, 'att', {'char','string'}}
                {A_bar, 'att', 'cell', 'nonempty'}
                {B_bar, 'att', 'cell', 'nonempty'}
                {dt, 'att', 'numeric', 'scalar'}
                });

    else
        inputArgsCheck({ ...
                {A_bar, 'att', 'cell', 'nonempty'}
                {B_bar, 'att', 'cell', 'nonempty'}
                {dt, 'att', 'numeric', 'scalar'}
                });
    end

end

end

function [name,A_bar,B_bar,dt,n_p,inputs,outputs] = ...
            aux_computeProperties(name,A_bar,B_bar,dt)

    % number of past time steps considered
    n_p = length(A_bar);

    % number of inputs, and outputs
    inputs = 1;
    if ~isscalar(B_bar{1,1})
        inputs = size(B_bar{1,1},2);
    end
    outputs = 1;
    if ~isscalar(A_bar{1,1})
        outputs = size(A_bar{1,1},1);
    end

end

% ------------------------------ END OF CODE ------------------------------
