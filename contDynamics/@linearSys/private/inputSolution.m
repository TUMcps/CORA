function [obj,options] = inputSolution(obj,options)
% inputSolution - computes the bloating due to the input 
%
% Syntax:
%    [obj] = inputSolution(obj,options)
%
% Inputs:
%    obj - linearSys object
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linearSys object
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       08-May-2007 
% Last update:   26-October-2007
%                07-October-2008
%                27-April-2009
%                22-July-2009
%                10-August-2010
%                15-June-2016
%                25-July-2016 (intervalhull replaced by interval)
%                01-November-2017 (consideration of constant input c)
%                16-November-2021 (disturbance set W)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% possible effect of disturbance: total input is then
%   V = obj.B*options.U + (options.W - center(options.W))
%   vTrans = obj.B*options.uTrans + center(options.W)
if isfield(options,'W')
    Wcenter = center(options.W);
    W = options.W + (-Wcenter);
else
    Wcenter = 0;
    W = 0;
end

% set of possible inputs
V = obj.B*options.U + W;

% do we have a time-varying input solution?
options.isRV = ~representsa_(V,'origin',eps);

% compute vTrans (including disturbance center and constant offset)
vTrans = obj.B*options.uTrans + Wcenter + obj.c;

Apower = obj.taylor.powers;
E = obj.taylor.error;
r = options.timeStep;
n = obj.dim;
factors = options.factor;

if options.isRV
    %init Vsum
    Vsum = r*V;
    Asum = r*eye(n);
    %compute higher order terms
    for i = 1:options.taylorTerms
        %compute sums
        Vsum = Vsum+Apower{i}*factors(i+1)*V;
        Asum = Asum+Apower{i}*factors(i+1);
    end
    
    %compute overall solution
    try
        inputSolV = Vsum+E*r*V;
    catch
        % for all set representations, which currently do not support the
        % (left-)multiplication with an interval matrix (in this case: E)
        inputSolV = Vsum+E*r*interval(V);
    end
    
else
    % only Asum, since V == origin (0)
    Asum = r*eye(n);
    %compute higher order terms
    for i = 1:options.taylorTerms
        %compute sum
        Asum = Asum+Apower{i}*factors(i+1);
    end
    inputSolV = zeros(n,1);
    
end

%compute solution due to constant input
eAtInt = Asum+E*r;
inputSolVtrans = eAtInt*zonotope(vTrans);

%compute additional uncertainty if origin is not contained in input set
if options.originContained
    inputCorr = zeros(n,1);
else
    %compute inputF
    obj = inputTie(obj,options);
    inputF = obj.taylor.inputF;
    inputCorr = inputF*zonotope(vTrans);
end


%write to object structure
obj.taylor.V = V;
obj.taylor.RV = inputSolV;
obj.taylor.Rtrans = inputSolVtrans;
obj.taylor.inputCorr = inputCorr;
obj.taylor.eAtInt = eAtInt;

% ------------------------------ END OF CODE ------------------------------
