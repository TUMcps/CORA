function [res,cert,scaling] = priv_zonotopeContainment_geneticMaximization(Z1, Z2, tol, maxEval, scalingToggle)
% priv_zonotopeContainment_geneticMaximization - Solves the zonotope containment problem by
%    checking whether the maximum value of the Z2-norm at one of the
%    vertices of Z1 exceeds 1+tol, using surrogate optimization,
%    see also [1].
%
% Syntax:
%    [res,cert,scaling] = priv_zonotopeContainment_geneticMaximization(Z1, Z2, tol, maxEval, scalingToggle)
%
% Inputs:
%    Z1 - zonotope object, inbody
%    Z2 - zonotope object, circumbody
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of Z
%       will be detected as lying in Z, which can be useful to counteract
%       errors originating from floating point errors.
%    maxEval - Number of maximal function evaluations.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below).
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, Z1 is
%           guaranteed to not be contained in Z2, whereas if res=false and
%           cert=false, nothing can be deduced (Z1 could still be
%           contained in Z2).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(Z2 - center(Z2)) + center(Z2) contains Z1.
%           For priv_zonotopeContainment_geneticMaximization, this is a lower bound.
%           Note that computing this scaling factor may significantly 
%           increase the runtime.
%
% Example: 
%    Z1 = zonotope([0.5 2 3 0;0.5 2 0 3]);
%    Z2 = zonotope([0 -1 1 0; 0 1 0 1]);
%    Z3 = Z2 + [3;0];
%    
%    % The function priv_zonotopeContainment_geneticMaximization will implicitly be called by
%    % contains
%    contains(Z1,Z2,'opt')
%    contains(Z1,Z3,'opt')
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
%    
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z3,[1,2],'r');
%
% References:
%    [1] A. Kulmburg, M. Althoff.: On the co-NP-Completeness of the
%        Zonotope Containment Problem, European Journal of Control 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/contains_

% Authors:       Adrian Kulmburg
% Written:       05-February-2024 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Retrieve the generator matrix of Z1 and determine its size
G = Z1.generators;
m = size(G, 2);

norm_Z2_nu = @(nu) -zonotopeNorm(Z2, G*nu' + Z1.center-Z2.center);

    function [state,options,optchanged] = customStoppingCriterion(options,state,flag)
        if state.FunEval > maxEval
            state.StopFlag = 'y';
        end
        optchanged = false;
    end

if scalingToggle
    persistent options_zonotope_contains___ga_scalingToggle_true
    if isempty(options_zonotope_contains___ga_scalingToggle_true)
        options_zonotope_contains___ga_scalingToggle_true = optimoptions("ga", ...
        'Display', 'iter', ... % Output on command line is iterative
        "PlotFcn", 'gaplotbestf',... % Plot value
        'OutputFcn', @customStoppingCriterion,...
        "MaxGenerations", inf);
    end
    options = options_zonotope_contains___ga_scalingToggle_true;
else
    persistent options_zonotope_contains___ga_scalingToggle_false
    if isempty(options_zonotope_contains___ga_scalingToggle_false)
        options_zonotope_contains___ga_scalingToggle_true = optimoptions("ga",...
        "Display", 'iter',... % Output on command line is iterative
        "PlotFcn", 'gaplotbestf',... % Plot value
        "FitnessLimit", -1-tol,...
        'OutputFcn', @customStoppingCriterion,...
        "MaxGenerations", inf);
    end
    options = options_zonotope_contains___ga_scalingToggle_false;
end

[~, scaling] = ga(norm_Z2_nu, m, [], [], [], [], -ones([m 1])', ones([m 1])', [], ones([m 1])', options);

scaling = abs(scaling);

if scaling > 1 + tol
    res = false;
    cert = true;
else
    res = true;
    cert = false;
end
    
end

% ------------------------------ END OF CODE ------------------------------
