classdef nnTanhLayer < nnSShapeLayer
% nnTanhLayer - class for tanh layers
%
% Syntax:
%    obj = nnTanhLayer(name)
%
% Inputs:
%    name - name of the layer, defaults to type
%
% Outputs:
%    obj - generated object
%
% References:
%    [1] Koller, L. "Co-Design for Training and Verifying Neural Networks",
%           Master's Thesis
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Tobias Ladner, Lukas Koller
% Written:       28-March-2022
% Last update:   25-May-2023 (LK, approx error: closed form for order=1)
% Last revision: 10-August-2022 (renamed)

% ------------------------------ BEGIN CODE -------------------------------

properties
    % found via obj.findRegionPolys
    reg_polys = [struct('l', -Inf,'u', -5,'p', [ 0.0000000000000000 -0.9999546021312975 ],'d', 0.0000453978687024),struct('l', -5,'u', -2.5000000000000000,'p', [ 0.0003632483926582 0.0077740610313991 0.0669424281432756 0.2908011544996281 0.6397668761140367 -0.4269522964773172 ],'d', 0.0000400131387163),struct('l', -2.5000000000000000,'u', -0.6250000000000000,'p', [ -0.0067465967733122 -0.0683314957855584 -0.2607942809217674 -0.3930263398561141 0.1220772321097194 1.1012160065134771 0.0237533346521373 ],'d', 0.0000714657357249),struct('l', -0.6250000000000000,'u', 1.2500000000000000,'p', [ -0.0055944532555784 -0.0272531518079237 0.1104113732706064 0.0164884302761415 -0.3299081428777030 -0.0024990843822790 0.9998499350506984 0.0000743310615162 ],'d', 0.0000967788202835),struct('l', 1.2500000000000000,'u', 5,'p', [ 0.0000512035784496 -0.0014556369685229 0.0178442098575226 -0.1227508298674645 0.5143045701658820 -1.3206066384526620 1.9385133540083290 -0.2654008527398747 ],'d', 0.0000712869080566),struct('l', 5,'u', Inf,'p', [ 0.0000000000000000 0.9999546021312975 ],'d', 0.0000453978687024)] 
end

methods
    % constructor
    function obj = nnTanhLayer(name)
        if nargin < 1
            name = [];
        end
        % call super class constructor
        obj@nnSShapeLayer(name)
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})
    % numeric
    function [r, obj] = evaluateNumeric(obj, input, options)
        r = tanh(input);
    end
end

methods
    function [coeffs, d] = computeApproxError(obj, l, u, coeffs)
        order = length(coeffs)-1;
        if order == 1
            % closed form: see [1]

            m = coeffs(1);
            t = coeffs(2);
            
            % compute extreme points of tanh - mx+t; there are two
            % solutions: xu and xl
            xu = atanh(sqrt(1 - m)); % point with max. upper error
            xl = -atanh(sqrt(1 - m)); % point with max. lower error

            % evaluate candidate extreme points within boundary
            xs = [l,xu,xl,u];
            xs = xs(l <= xs & xs <= u);
            ys = obj.f(xs);
            
            % compute approximation error at candidates
            dBounds = ys - (m*xs+t);
            
            % compute approximation error
            du = max(dBounds);
            dl = min(dBounds);
            dc = (du+dl)/2;
            d = du-dc;

            % shift coeffs by center
            coeffs = [m, t + dc];

            computeApproxError@nnActivationLayer(obj,l,u,coeffs);
        else
            % compute in super class
            [coeffs,d] = computeApproxError@nnSShapeLayer(obj,l,u,coeffs);
        end
    end
end

methods(Access=protected)
    function [xs,dxsdm] = computeExtremePointsBatch(obj, m, options)
        % compute extreme points of tanh - mx+t; there are two
        % solutions: xu and xl
        m = max(min(1,m),eps('like',m));
        xu = atanh(sqrt(1 - m)); % point with max. upper error
        xl = -xu; % point with max. lower error
        % list of extreme points
        xs = cat(3,xl,xu);
        % compute derivate wrt. slope m; needed for backpropagation
        dxu = -1./max(2*m.*sqrt(1 - m),eps('like',m));
        dxl = -dxu;
        dxsdm = cat(3,dxl,dxu);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
