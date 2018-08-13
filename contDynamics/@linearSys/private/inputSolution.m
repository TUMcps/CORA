function [obj] = inputSolution(obj,options)
% inputSolution - computes the bloating due to the input 
%
% Syntax:  
%    [obj] = inputSolution(obj,options)
%
% Inputs:
%    obj - linearSys object
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

% Author: Matthias Althoff
% Written: 08-May-2007 
% Last update:  26-October-2007
%               07-October-2008
%               27-April-2009
%               22-July-2009
%               10-August-2010
%               15-June-2016
%               25-July-2016 (intervalhull replaced by interval)
%               01-November-2017 (consideration of constant input c)
% Last revision:---

%------------- BEGIN CODE --------------

% set of possible inputs
V = obj.B*options.U;

% compute vTrans 
vTrans = obj.B*options.uTrans;

% consider constant input
if ~isempty(obj.c)
    vTrans = vTrans + obj.c;
end

A = obj.A;
Apower = obj.taylor.powers;
E = obj.taylor.error;
taylorTerms = options.taylorTerms;
r = options.timeStep;
dim = length(A);


%init Vsum
Vsum = r*V;
Asum = r*eye(dim);
%compute higher order terms
for i = 1:taylorTerms
    %compute sums
    Vsum = Vsum+Apower{i}*r^(i+1)/factorial(i+1)*V;
    Asum = Asum+Apower{i}*r^(i+1)/factorial(i+1);
end

%compute overall solution
inputSolV = Vsum+E*r*V;

%compute solution due to constant input
eAtInt = Asum+E*r;
inputSolVtrans = eAtInt*zonotope(vTrans);

%compute additional uncertainty if origin is not contained in input set
if options.originContained
    inputCorr = zeros(dim,1);
else
    %compute inputF
    [obj] = inputTie(obj,options);
    inputF = obj.taylor.inputF;
    inputCorr = inputF*zonotope(vTrans);
end


%write to object structure
obj.taylor.V = V;
obj.taylor.RV = inputSolV;
obj.taylor.Rtrans = inputSolVtrans;
obj.taylor.inputCorr = inputCorr;
obj.taylor.eAtInt = eAtInt;

%------------- END OF CODE --------------