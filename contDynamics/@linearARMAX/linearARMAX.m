classdef linearARMAX < contDynamics
% linearARMAX - object constructor for linear discrete-time ARMAX systems
%
% Syntax:
%    obj = linearARMAX(A_bar,B_bar,dt)
%    obj = linearARMAX(name,A_bar,B_bar,dt)
%
% Description:
%    Generates a discrete-time linear ARMAX object according to the 
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
%    obj - generated linearARMAX object
%
% Example:
%    dt = 0.1;
%    A_bar = {[-0.4 0.6; 0.6 -0.4];[0.1 0; 0.2 -0.5]};
%    B_bar = {[0; 0];[0.3; -0.7];[0.1; 0]};
%    sys = linearARMAX(A_bar,B_bar,dt)
%
% Reference:
%   [1] L. Luetzow, M. Althoff. "Reachability Analysis of ARMAX Models", in
%       Proc. of the 62th International Conference on Decision and Control,
%       2023.

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
end

methods

    % class constructor
    function obj = linearARMAX(varargin)

        % not enough or too many input arguments
        if nargin < 3
            throw(CORAerror('CORA:notEnoughInputArgs',3));
        elseif nargin > 4
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end

        % parse name, output parameters, input parameters
        if ischar(varargin{1})
            name = varargin{1};
            varargin = varargin(2:end);
        else
            name = 'linearARMAX'; % default name
        end

        A_bar = varargin{1};
        B_bar = varargin{2};
        dt = varargin{3};

        % number of past time steps considered
        p = length(A_bar);

        % number of inputs, and outputs
        inputs = 1;
        outputs = 1;

        if ~isscalar(B_bar{1,1})
            inputs = size(B_bar{1,1},2);
        end
        if ~isscalar(A_bar{1,1})
            outputs = size(A_bar{1,1},1);
        end
        inputArgsCheck({ ...
            {A_bar, 'att', 'cell', 'nonempty'}
            {B_bar, 'att', 'cell', 'nonempty'}
            {dt, 'att', 'numeric', 'scalar'}
            });

        % instantiate parent class
        obj@contDynamics(name,p,inputs,outputs);

        % assign object properties
        obj.A_bar = A_bar;
        obj.B_bar = B_bar;
        obj.dt = dt;
        obj.tvp = false;
    end

    % invoke function observe so that the superclass can access private
    % functions
    function [R, tcomp] = observe(obj,params,options)
        [R, tcomp] = observe@contDynamics(obj,params,options);
    end

    function sys = setTVP(sys)
        % set tvp attribute to true and compute the parameters required for 
        % the calculation of the time-varyig parameters A_tilde and B_tilde

        sys.tvp = true;
        if ~isfield(sys, "conv_tvp")
            p = sys.dim;
            dim_y = sys.nrOfOutputs;

            % build augmented matrices
            sys.conv_tvp.A_ext = [zeros(dim_y*(p-1), dim_y) eye(dim_y * (p-1)); zeros(dim_y, dim_y*p)];
            sys.conv_tvp.B_ext = cell(p+1,1);
            sys.conv_tvp.E = [zeros(dim_y, dim_y * (p-1)) eye(dim_y)];
            for i = 0:p
                if i < p
                    col_i = i*dim_y+1:(i+1)*dim_y;
                    sys.conv_tvp.A_ext(end-dim_y+1:end,col_i) = sys.A_bar{p-i};
                end
                sys.conv_tvp.B_ext{i+1} = sys.conv_tvp.E'*sys.B_bar{i+1};
            end
        end
    end

    function sys = negateTVP(sys)
        % set tvp attribute to false

        sys.tvp = false;
    end

    function sys = computeTVP(sys, k, iVec)
        % compute time-varying parameters A_tilde and B_tilde for the time 
        % point k and the indezes in iVec

        p = sys.dim;
        k_plus = k+p-1;
        dim_u = sys.nrOfInputs;
        dim_y = sys.nrOfOutputs;
        if ~sys.tvp
            sys.setTVP;
        end

        if isfield(sys, "A_tilde")
            k_last = length(sys.A_tilde,1);
        else
            k_last = 0;
            sys.A_tilde = cell(k+1,1);
            sys.B_tilde = cell(max(iVec)+1,k+1);
        end

        if max(iVec) > k_plus || min(iVec) < 0
            error("Invalid. Index i must be in 0<=i<=k_plus.")
        end

        if k<p
            error("Invalid. Time step k must be in k>=p.")
        elseif k == 1
            sys.A_tilde{k+1} = sys.conv_tvp.A_ext;
            for i = iVec
                sys.B_tilde{i+1,k+1} = sys.conv_tvp.B_ext{i+1};
            end
        else 
            if k_last == 0 
                sys.A_tilde{k+1} = sys.conv_tvp.A_ext^(k); % [1, Eq. (9a)] 
            else
                sys.A_tilde{k+1} = sys.A_tilde{k_last} * ...
                    sys.conv_tvp.A_ext^(k-k_last); % [1, Eq. (14a)] 
            end

            % compute parameters from scratch (see [1, Eq. (9b)])
            for i = iVec
                j_min = max(0,i-p);
                j_max = min(k-1, i);
                sys.B_tilde{i+1,k+1} = zeros(p*dim_y, dim_u);
                Aj = sys.conv_tvp.A_ext^j_min;
                for j=j_min:j_max
                    if j > j_min 
                        Aj = Aj * sys.conv_tvp.A_ext;
                    end
                    sys.B_tilde{i+1,k+1} = sys.B_tilde{i+1,k+1} + Aj * sys.conv_tvp.B_ext{i-j+1};
                end
            end
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
