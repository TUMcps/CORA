function x_ = reset(trans,x,varargin)
% reset - resets the continuous state according the reset function
%
% Syntax:
%    x_ = reset(trans,x)
%    x_ = reset(trans,x,u)
%
% Inputs:
%    trans - transition object
%    x - state value before the reset (vector, set, or cell thereof)
%    u - (optional) input at the time of reset (vector, set, or cell thereof)
%
% Outputs:
%    x_ - state value after the reset
%
% Example:
%    trans = transition(conHyperplane([1 0],0),...
%            struct('A',eye(2),'c',[1;-1]),1);
%    x = [2;1];
%    x_ = reset(trans,x);
%
% References:
%    [1] N. Kochdumper and et al. "Reachability Analysis for Hybrid Systems
%        with Nonlinear Guard Sets", HSCC 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: transition

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       04-May-2007 
% Last update:   07-October-2008
%                10-December-2021 (NK, added nonlinear reset functions)
%                11-February-2022 (MP, added input dependence)
%                11-May-2023 (MW, recursive call for cells, bug fixes)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check whether multiple points/sets given
if iscell(x)
    % default: no inputs given
    u = cell(length(x),1);
    % check if inputs are given
    if nargin > 2
        u = varargin{1};
    end

    % recursive call
    x_ = cell(length(x),1);
    for i=1:length(x_)
        x_{i} = reset(trans,x{i},u{i});
    end

else

    % default: no inputs given
    u = [];
    % check if inputs were given
    if nargin > 2
        u = varargin{1};
    end

    % get correct orientation of x
    if size(x,2) > 1
        x = x';
    end

    % check whether reset function is linear or nonlinear
    if isfield(trans.reset,'A')
        % linear reset
        if ~isfield(trans.reset,'B') || (~any(any(trans.reset.B)) && isempty(u))
            % without inputs
            x_ = trans.reset.A*x + trans.reset.c;
        else
            % with inputs
            x_ = trans.reset.A*x + trans.reset.B*u + trans.reset.c;
        end

    else
        % nonlinear reset
        if trans.reset.inputDim > 0
            % compute extended state x' = [x;u]
            x = cartProd_(x,u,'exact');
        end
        if isnumeric(x)
            x_ = trans.reset.f(x);
        else
            x_ = aux_nonlinearReset(trans,x);
        end
    end

end

end


% Auxiliary functions -----------------------------------------------------

function X_ = aux_nonlinearReset(trans,Z)
% compute the resulting set for a nonlinear reset function using the
% approach described in Sec. 4.4 in [1]

    % dimension of state and (global) inputs
    n = trans.reset.stateDim;
    m = trans.reset.inputDim;

    % Taylor series expansion point
    % (length: number of states + number of global inputs)
    p = center(Z);
    
    % interval enclosure of set
    I = interval(Z);
    
    % evaluate reset function at expansion point
    f = trans.reset.f(p);

    % evaluate Jacobian at expansion point
    J = trans.reset.J(p);
    
    % evaluate Hessians at expansion point
    Q = repmat({0},[n,1]);
    for i = 1:n
        funHan = trans.reset.Q{i};
        if ~isnumeric(funHan)
            % if Hessian is numeric, it is zero
            Q{i} = funHan(p);
        end
    end
    
    % evaluate third-order tensors over entire set
    T = repmat({0},[n,n+m]);
    for i = 1:n
        for j = 1:n+m
            funHan = trans.reset.T{i,j};
            if ~isnumeric(funHan)
                % if third-order tensor is numeric, it is zero
                T{i,j} = funHan(I);
            end
        end
    end

    % compute Lagrange remainder
    I = I - p;
    rem = interval(zeros(n,1));
    
    for i = 1:n
        for j = 1:n+m
            if ~isnumeric(T{i,j})
                % skip entries that do not contribute to remainder sum
                rem(i) = rem(i) + I(j) * transpose(I) * T{i,j} * I;
            end
        end
    end
    
    % include factor and convert remainder to zonotope for quadMap below
    rem = zonotope(1/6*rem);

    % compute over-approximation using a Taylor series of order 2
    X_ = f + J * (Z + (-p)) + 0.5*quadMap((Z + (-p)),Q) + rem;
end

% ------------------------------ END OF CODE ------------------------------
