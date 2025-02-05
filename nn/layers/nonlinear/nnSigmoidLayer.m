classdef nnSigmoidLayer < nnSShapeLayer
% nnSigmoidLayer - class for Sigmoid layers
%
% Syntax:
%    obj = nnSigmoidLayer(name)
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
    reg_polys = [ struct('l', -Inf,'u', -10,'p', [ 0.0000000000000000 0.0000226989343512 ],'d', 0.0000226989343512),struct('l', -10,'u', -5,'p', [ 0.0000056757561353 0.0002429394072305 0.0041839017589345 0.0363501443122013 0.1599417190270924 0.2865242403387572 ],'d', 0.0000112339935969),struct('l', -5,'u', -1.2500000000000000,'p', [ -0.0000527077872915 -0.0010676796216488 -0.0081498212788008 -0.0245641462409732 0.0152596540138349 0.2753040016287754 0.5118757437687305 ],'d', 0.0000267218430910),struct('l', -1.2500000000000000,'u', 2.5000000000000000,'p', [ 0.0000281300436054 -0.0001624932291702 -0.0000940034934591 0.0020483035259502 0.0001012481844055 -0.0208261173225880 -0.0000360364456032 0.2499994775884087 0.5000017612516241 ],'d', 0.0000093721843051),struct('l', 2.5000000000000000,'u', 10,'p', [ 0.0000000076826740 -0.0000001836711225 -0.0000031954035029 0.0001818834877788 -0.0031382230953991 0.0290271690503697 -0.1566526667804534 0.4720559527908412 0.3752299739623118 ],'d', 0.0000130525253081),struct('l', 10,'u', Inf,'p', [ 0.0000000000000000 0.9999773010656487 ],'d', 0.0000226989343513)]
end

methods
    % constructor
    function obj = nnSigmoidLayer(name)
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
    function r = evaluateNumeric(obj, input, options)
        % use tanh for numeric stability
        r = tanh(input/2) / 2 + 0.5;
    end
end

methods
    function [coeffs, d] = computeApproxError(obj, l, u, coeffs)
        order = length(coeffs)-1;
        if order == 1
            % closed form: see [1]

            m = coeffs(1);
            t = coeffs(2);

            % compute extreme points of sigmoid - mx+t; there are two
            % solutions: xu and xl
            xu = 2*atanh(sqrt(1 - 4*m)); % point with max. upper error
            xl = -2*atanh(sqrt(1 - 4*m)); % point with max. lower error

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
        else
            % compute in super class
            [coeffs,d] = computeApproxError@nnSShapeLayer(obj,l,u,coeffs);
        end
    end
end

methods(Access=protected)
    function [xs,dxsdm] = computeExtremePointsBatch(obj, m, options)
        % compute extreme points of sigmoid - mx+t; there are two
        % solutions: xu and xl
        m = min(1/4,m);
        xu = 2*atanh(sqrt(1 - 4*m)); % point with max. upper error
        xl = -xu; % point with max. lower error
        % list of extreme points
        xs = cat(3,xl,xu);
        % compute derivate wrt. slope m; needed for backpropagation
        dxu = -1./max(m.*sqrt(1 - 4*m),eps('like',m));
        dxl = -dxu;
        dxsdm = cat(3,dxl,dxu);
    end
end

end

% ------------------------------ END OF CODE ------------------------------
