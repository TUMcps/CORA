function coeffs = findRegionPolys(obj,tol,order,l_max,u_max,pStart,dStart,pEnd,dEnd)
% findRegionPolys - divides the given domain of the nonlinear function into 
%     regions and finds polynomials such that the approximation error is 
%     at most tol
%
% Syntax:
%    coeffs = findRegionPolys(tol,l_max,u_max,coeffStart,dStart,coeffEnd,dEnd)
%
% Inputs:
%    obj - nnActivationLayer
%    tol - tolerance
%    order - order of polynomials
%    l_max,u_max - bounds for domain
%    pStart - polynom for [-inf,l_max]
%    dStart - error of pStart
%    pEnd - polynom for [u_max,inf]
%    dEnd - error of pEnd
%
% Outputs:
%    coeffs - struct with
%        - .l - lower bound of active region
%        - .u - upper bound of active region
%        - .p - polynomial coefficients
%        - .d - max error
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner
% Written:       26-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin < 9
    throw(CORAerror('CORA:notEnoughInputArgs', 9))
elseif nargin > 9
    throw(CORAerror('CORA:tooManyInputArgs', 9))
end

inputArgsCheck({ ...
    {obj, 'att', 'nnActivationLayer'}, ...
    {tol, 'att', 'numeric', 'scalar'}, ...
    {order, 'att', 'numeric', 'scalar'}, ...
    {l_max, 'att', 'numeric', 'scalar'}, ...
    {u_max, 'att', 'numeric', 'scalar'}, ...
    {pStart, 'att', 'numeric', 'row'}, ...
    {dStart, 'att', 'numeric', {'scalar','positive'}}, ...
    {pEnd, 'att', 'numeric', 'row'}, ...
    {dEnd, 'att', 'numeric', {'scalar','positive'}}, ...
})

if dStart > tol
    warning('CORA: nnActivationLayer/findRegionPolys: dStart > tol.')
end
if dEnd > tol
    warning('CORA: nnActivationLayer/findRegionPolys: dEnd > tol.')
end

% remove reg_polys if present
if isprop(obj, 'reg_polys')
    warning("CORA: Temporarily removing current region polynomials.")
    reg_polys = obj.reg_polys;
    obj.reg_polys = [];
end

% init struct
coeffs = struct;

% start polynomial
coeffs(1).l = -inf;
coeffs(1).u = l_max;
coeffs(1).p = pStart;
coeffs(1).d = dStart;

% iteratively find approx polynomial
l = l_max;
u = u_max;
while l < u
    fprintf("Progress: %.2f%%\n", ((u_max-l_max)-(u-l))/(u_max-l_max) * 100)

    % left polynomial ---
    u_i = u;
    d = Inf;
    p = [0 0];

    while d > tol
        if ~isinf(d)
            % shorten region
            u_i = l + (u_i-l) * 0.5;

            if (u_i-l) < tol
                throw(CORAerror("CORA:notConverged",sprintf('Unable to find polynomial at [%.6f,%.6f]. Reached x=%.6f. y',l,u,u_i)))
            end
        end

        % find polynomial + error
        [p,d] = obj.computeApproxPoly(l,u_i,order,"regression");
    end

    % working polynomial found, trying to reduce the order
    for o=order-1:-1:0
        % find polynomial + error
        [p_o,d_o] = obj.computeApproxPoly(l,u_i,o,"regression");

        if d_o < tol
            % update with lower-order polynomail
            p = p_o;
            d = d_o;
        else
            % no lower order polynomial will recover d to be below tol
            break;
        end


    end

    % next polynomial
    coeffs(end+1).l = l;
    coeffs(end).u = u_i;
    coeffs(end).p = p;
    coeffs(end).d = d;

    % move to next region
    l = u_i;

    if l == u
        % check if entire domain is covered
        break
    end
    
    fprintf("Progress: %.2f%%\n", ((u_max-l_max)-(u-l))/(u_max-l_max) * 100)

    % right polynomial ---
    l_i = l;
    d = Inf;
    p = [0 0];

    while d > tol
        if ~isinf(d)
            % half region
            l_i = u - (u-l_i)/2;

            if (u-l_i) < tol
                throw(CORAerror("CORA:notConverged",sprintf('Unable to find polynomial at [%.6f,%.6f]. Reached x=%.6f. y',l,u,l_i)))
            end
        end

        % find polynomial + error
        p = obj.computeApproxPoly(l_i,u,order,"regression");
        [p,d] = obj.computeApproxError(l_i,u,p);
    end

    % next polynomial
    coeffs(end+1).l = l_i;
    coeffs(end).u = u;
    coeffs(end).p = p;
    coeffs(end).d = d;

    % move to next region
    u = l_i;
end
fprintf("Progress: %.2f%%\n", ((u_max-l_max)-(u-l))/(u_max-l_max) * 100)

% end polynomial
coeffs(end+1).l = u_max;
coeffs(end).u = inf;
coeffs(end).p = pEnd;
coeffs(end).d = dEnd;

% sort coeffs
[~,idx] = sort([coeffs.l]);
coeffs = coeffs(idx);

fprintf("Required %d polynomials.\n", length(coeffs));

% restore reg_polys
% remove reg_polys if present
if isprop(obj, 'reg_polys')
    obj.reg_polys = reg_polys;
end

end

% ------------------------------ END OF CODE ------------------------------
