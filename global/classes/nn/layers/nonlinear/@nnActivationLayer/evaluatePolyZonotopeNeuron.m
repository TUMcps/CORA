function [c, G, GI, d] = evaluatePolyZonotopeNeuron(obj, c, G, GI, E, E_out, order, ind, ind_, evParams)
% evaluatePolyZonotopeNeuron - evaluates a single activation neuron on a 
%    polyZonotope
%
% Syntax:
%    [c, G, GI, d] = evaluatePolyZonotopeNeuron(obj, c, G, GI, E, Es, id, id_, ind, ind_, evParams)
%
% Inputs:
%    c, G, GI, E, id, id_, ind, ind_ - parameters of input polyZonotope
%    E_out - exponent matrix of the output polyZonotope
%    order - polynomial degree used for approximation
%    evParams - parameter for NN evaluation
%
% Outputs:
%    updated [c, G, GI, d]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnActivationLayer/evaluatePolyZonotope

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   05-April-2022
%                13-May-2022 (reuse bounds)
%                24-June-2022 (performance optimizations)
%                05-December-2022 (readability through aux functions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enclose the activation function with a polynomial zonotope
% by fitting a polynomial function

% compute lower and upper bound of domain
[l, u] = aux_computeBounds(obj, c, G, GI, E, ind, ind_, evParams);

% compute approximating polynomial 
[coeffs, d] = computeApproxPoly(obj, l, u, order, evParams.poly_method);

% evaluate the polynomial approximation on the polynomial zonotope
[c, G, GI] = aux_evaluatePolyZonotopeOnPolynomial(obj, c, G, GI, E_out, coeffs, order);

% approx error is added after all neurons are processed
end


% Auxiliary functions -----------------------------------------------------

function [l, u] = aux_computeBounds(obj, c, G, GI, E, ind, ind_, evParams)
    % compute lower and upper bound of domain

    if ~evParams.reuse_bounds || isnan(obj.l(evParams.j)) || isnan(obj.u(evParams.j))
        approx = evParams.bound_approx;
        if isa(approx, 'logical')
            [l, u] = nnHelper.compBoundsPolyZono(c, G, GI, E, ind, ind_, approx);
        elseif strcmp(approx, "sample")
            % not used anymore 
            pZ = polyZonotope(c, G, GI, E);
            points = [pZ.randPoint(5000), pZ.randPoint(500, 'extreme')];
            tol = 0.00001;
            int = interval(min(points), max(points));
            int = int + int * tol;
            l = infimum(int);
            u = supremum(int);
        else
            throw(CORAerror("CORA:wrongFieldValue", ...
                "evParams.bound_approx", {true, false, 'sample'}))
        end
    
        if evParams.reuse_bounds
            obj.l(evParams.j) = l;
            obj.u(evParams.j) = u;
        end
    else
        % get quick estimate
        [l, u] = nnHelper.compBoundsPolyZono(c, G, GI, E, ind, ind_, true);
        
        % determine tighter bounds
        obj.l(evParams.j) = max(obj.l(evParams.j), l);
        obj.u(evParams.j) = min(obj.u(evParams.j), u);
    
        % reuse bounds
        l = obj.l(evParams.j);
        u = obj.u(evParams.j);
    end

    if l > u 
        if withinTol(l, u)
            % fix numerical instability
            l_ = l;
            l = min(l, u);
            u = max(l_, u);

            obj.l(evParams.j) = l;
            obj.u(evParams.j) = u;
        else
            warnText = 'Layer %i, neuron %i: Lower bound is larger than upper bound.';
            if evParams.reuse_bounds
                warnText = [warnText, ' Might be due to reused bounds.'];
            end
            throw(CORAerror('CORA:specialError', sprintf( ...
                warnText, evParams.i, evParams.j ...
            )))
        end
    end
end

function [c, G, GI] = aux_evaluatePolyZonotopeOnPolynomial(obj, c, G, GI, E_, coeffs, order)
    % evaluate the polynomial approximation on the polynomial zonotope

    % higher-order might not be required if activation function can be exploited
    % e.g. fully in linear segment for ReLU
    order = length(coeffs)-1;
    
    % store information to directly compact generators
    [G_start, G_end, G_ext_start, G_ext_end] = nnHelper.getOrderIndicesG(G, order);
    [GI_start, GI_end, GI_ext_start, GI_ext_end] = nnHelper.getOrderIndicesGI(GI, G, order);

    % init
    c_out = zeros(1, order);
    G_out = zeros(1, G_end(end));
    G_ext_out = zeros(1, G_ext_end(end));
    GI_out = zeros(1, GI_end(end));
    GI_ext_out = zeros(1, GI_ext_end(end));
    
    c_out(1) = c;
    G_ext_out(G_start(1):G_end(1)) = G;
    GI_ext_out(GI_start(1):GI_end(1)) = GI;
    
    for i = 2:order
        % Note that e.g., G2 * G3 = G5
        i1 = floor(i/2);
        i2 = ceil(i/2);
    
        % previous values
        ci1 = c_out(i1);
        Gi1 = G_ext_out(G_ext_start(i1):G_ext_end(i1));
        GIi1 = GI_ext_out(GI_ext_start(i1):GI_ext_end(i1));
        Ei1 = E_(:, G_ext_start(i1):G_ext_end(i1));
    
        ci2 = c_out(i2);
        Gi2 = G_ext_out(G_ext_start(i2):G_ext_end(i2));
        GIi2 = GI_ext_out(GI_ext_start(i2):GI_ext_end(i2));
        Ei2 = E_(:, G_ext_start(i2):G_ext_end(i2));
    
        % calculate i
        [ci, Gi_ext, GIi_ext] = nnHelper.calcSquared(ci1, Gi1, GIi1, Ei1, ci2, Gi2, GIi2, Ei2, i1 == i2);
    
        % store results
        c_out(i) = ci;
        G_ext_out(G_ext_start(i):G_ext_end(i)) = Gi_ext;
        GI_ext_out(GI_ext_start(i):GI_ext_end(i)) = GIi_ext;
    end
    
    % update weights with coefficients
    for i = 1:order
        coeff_i = coeffs(end-i);
        c_out(i) = c_out(i) * coeff_i;
        G_ext_out(G_ext_start(i):G_ext_end(i)) = G_ext_out(G_ext_start(i):G_ext_end(i)) * coeff_i;
        GI_ext_out(GI_ext_start(i):GI_ext_end(i)) = GI_ext_out(GI_ext_start(i):GI_ext_end(i)) * coeff_i;
    end
    
    % update generators with same exponent from back to front:
    % i1 = floor(i/2); i2 = ceil(i/2)
    % G{i} = [Gi1, Gi2, Gi]
    % → G{i1} += Gi1
    % → G{i2} += Gi2
    % and add Gi to result
    % analogous for GI{i} = [GIi1, GIi2, GIi]
    
    for i = order:-1:2
        % Note that e.g., G2 * G3 = G5
        i1 = floor(i/2);
        i2 = ceil(i/2);
    
        % get generator lenghts
        i_len = G_end(i) - G_start(i) + 1;
        i1_len = G_ext_end(i1) - G_ext_start(i1) + 1;
        i2_len = G_ext_end(i2) - G_ext_start(i2) + 1;
    
        % extract lower order generators
        Gi1 = G_ext_out(1:end, G_ext_start(i)-1+(1:i1_len));
        Gi2 = G_ext_out(1:end, G_ext_start(i)+i1_len-1+(1:i2_len));
        Gi = G_ext_out(1:end, G_ext_start(i)+i1_len+i2_len-1+(1:i_len));
    
        % update generators
        G_ext_out(G_ext_start(i1):G_ext_end(i1)) = G_ext_out(G_ext_start(i1):G_ext_end(i1)) + Gi1;
        G_ext_out(G_ext_start(i2):G_ext_end(i2)) = G_ext_out(G_ext_start(i2):G_ext_end(i2)) + Gi2;
        G_out(G_start(i):G_end(i)) = Gi;
    
        % analogous for GI if present
        if GI_ext_end(i) > 0
            % get generator lenghts
            i_len = GI_end(i) - GI_start(i) + 1;
            i1_len = GI_ext_end(i1) - GI_ext_start(i1) + 1;
            i2_len = GI_ext_end(i2) - GI_ext_start(i2) + 1;
    
            % extract lower order generators
            GIi1 = GI_ext_out(1:end, GI_ext_start(i)-1+(1:i1_len));
            GIi2 = GI_ext_out(1:end, GI_ext_start(i)+i1_len-1+(1:i2_len));
            GIi = GI_ext_out(1:end, GI_ext_start(i)+i1_len+i2_len-1+(1:i_len));
    
            % update lower order generators
            GI_ext_out(GI_ext_start(i1):GI_ext_end(i1)) = GI_ext_out(GI_ext_start(i1):GI_ext_end(i1)) + GIi1;
            GI_ext_out(GI_ext_start(i2):GI_ext_end(i2)) = GI_ext_out(GI_ext_start(i2):GI_ext_end(i2)) + GIi2;
            GI_out(GI_start(i):GI_end(i)) = GIi;
        end
    end
    
    % extract linear generators
    G_out(G_start(1):G_end(1)) = G_ext_out(G_ext_start(1):G_ext_end(1));
    GI_out(GI_start(1):GI_end(1)) = GI_ext_out(GI_ext_start(1):GI_ext_end(1));
    
    c = sum(c_out) + coeffs(end);
    G = G_out;
    GI = GI_out;

end


% ------------------------------ END OF CODE ------------------------------
