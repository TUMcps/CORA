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
    reg_polys = [ ...
     struct( ...
        'l',             -Inf, ...
        'u',              -10, ...
        'p', [                0, 2.269893435122294e-05 ], ...
        'd', 2.269893435122294e-05 ...
    ),struct( ...
        'l',              -10, ...
        'u',               -5, ...
        'p', [ 5.675756135273098e-06, 0.0002429394072304745, 0.004183901758934464, 0.03635014431220126, 0.1599417190270924, 0.2865242403387572 ], ...
        'd', 1.123399359691206e-05 ...
    ),struct( ...
        'l',               -5, ...
        'u',            -1.25, ...
        'p', [ -5.270778729152262e-05, -0.001067679621648762, -0.008149821278800829, -0.02456414624097322, 0.01525965401383487, 0.2753040016287754, 0.5118757437687305 ], ...
        'd', 2.672184309100989e-05 ...
    ),struct( ...
        'l',            -1.25, ...
        'u',              2.5, ...
        'p', [ 2.813004360543072e-05, -0.0001624932291701717, -9.400349345914858e-05, 0.002048303525950239, 0.0001012481844055254, -0.02082611732258803, -3.603644560320647e-05, 0.2499994775884087, 0.5000017612516241 ], ...
        'd', 9.372184305144976e-06 ...
    ),struct( ...
        'l',              2.5, ...
        'u',               10, ...
        'p', [ 7.682673983224953e-09, -1.836711225209388e-07, -3.195403502850288e-06, 0.0001818834877788483, -0.003138223095399062, 0.02902716905036973, -0.1566526667804534, 0.4720559527908412, 0.3752299739623118 ], ...
        'd', 1.305252530809456e-05 ...
    ),struct( ...
        'l',               10, ...
        'u',              Inf, ...
        'p', [                0, 0.9999773010656487 ], ...
        'd', 2.26989343512507e-05 ...
    ); ...
    ]
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
