function y = reset(trans,x,varargin)
% reset - resets the continuous state according the reset function
%
% Syntax:  
%    y = reset(trans,x)
%    y = reset(trans,x,u)
%
% Inputs:
%    trans - transition object
%    x - state value before the reset
%    u - (optional) input at the time of reset
%
% Outputs:
%    y - state value after the reset
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

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      04-May-2007 
% Last update:  07-October-2008
%               10-December-2021 (NK, added nonlinear reset functions)
%               11-February-2022 (MP, added input dependence)
% Last revision:---

%------------- BEGIN CODE --------------

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

    % loop over all points
    if ~iscell(x)

        % check whether reset function is linear or nonlinear
        if isfield(trans.reset,'A')
            % linear reset

            if ~isfield(trans.reset,'B') || (~any(any(trans.reset.B)) && isempty(u))
                % without inputs
                y = trans.reset.A*x + trans.reset.c;
            else
                % with inputs
                y = trans.reset.A*x + trans.reset.B*u + trans.reset.c;
            end

        else
            % nonlinear reset

            if trans.reset.inputDim > 0
                % compute extended state x' = [x;u]
                x = cartProd(x,u);
            end
            if isnumeric(x)
                y = trans.reset.f(x);
            else
                y = nonlinearReset(trans,x);
            end
        end

    else
        y = cell(length(x));
        for i=1:length(x)
            if ~isempty(x{i})
                if isfield(trans.reset,'A')
                    % linear reset
                    if isfield(trans.reset,'B')
                        % with input
                        y{i} = trans.reset.A*x{i} + trans.reset*B + trans.reset.c;
                    else
                        % without input
                        y{i} = trans.reset.A*x{i} + trans.reset.c;
                    end
                
                else
                    % nonlinear reset
                    x{i} = cartProd(x{i},u{i});
                    if isnumeric(y{i})
                        y{i} = trans.reset.f(x{i});
                    else
                        y{i} = nonlinearReset(trans,x{i});
                    end
                end

            end
        end
    end

end


% Auxiliary Functions -----------------------------------------------------

function Y = nonlinearReset(trans,X)
% compute the resulting set for a nonlinear reset function using the
% approach described in Sec. 4.4 in [1]

    % dimension of state and (global) inputs
    n = trans.reset.stateDim;
    m = trans.reset.inputDim;

    % Taylor series expansion point
    % (length: number of states + number of global inputs)
    p = center(X);
    
    % interval enclosure of set
    I = interval(X);
    
    % evaluate reset function at expansion point
    f = trans.reset.f(p);

    % evaluate Jacobian at expansion point
    J = trans.reset.J(p);
    
    % evaluate Hessians at expansion point
    Q = cell(n,1);
    for i = 1:n
        funHan = trans.reset.Q{i};
        if ~isempty(funHan)
            Q{i} = funHan(p);
        end
    end
    
    % evaluate third-order tensors over entire set
    T = cell(n,n+m);
    for i = 1:n
        for j = 1:n+m
            funHan = trans.reset.T{i,j};
            if ~isempty(funHan)
                T{i,j} = funHan(I);
            end
        end
    end

    % compute Lagrange remainder
    I = I - p;
    rem = interval(zeros(n,1));
    
    for i = 1:n
        for j = 1:n+m
            if ~isempty(T{i,j})
                rem(i) = rem(i) + I(j) * transpose(I) * T{i,j} * I;
            end
        end
    end
    
    % include factor and convert remainder to zonotope for quadMap below
    rem = zonotope(1/6*rem);

    % compute over-approximation using a Taylor series of order 2
    Y = f + J * (X + (-p)) + 0.5*quadMap((X + (-p)),Q) + rem;
end

%------------- END OF CODE --------------