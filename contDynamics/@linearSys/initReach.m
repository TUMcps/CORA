function [Rfirst,options] = initReach(linsys,Rinit,params,options)
% initReach - computes the reachable continuous set for the first time step
%
% Syntax:
%    [Rfirst,options] = initReach(linsys,Rinit,params,options)
%
% Inputs:
%    linsys - linearSys object
%    Rinit - initial reachable set
%    params - model parameters
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rfirst - first reachable set 
%    options - options for the computation of the reachable set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       07-May-2007 
% Last update:   03-January-2008
%                04-May-2009
%                29-June-2009
%                08-August-2011
%                25-July-2016 (intervalhull replaced by interval)
%                06-April-2017
%                28-October-2017
%                07-November-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute exponential matrix
linsys = aux_exponential(linsys,options);
% compute time interval error (tie)
linsys = aux_tie(linsys,options);
% compute reachable set due to input
linsys = aux_inputSolution(linsys,params,options);
%change the time step
linsys.taylor.timeStep=options.timeStep;

%compute reachable set of first time interval
eAt=expm(linsys.A*options.timeStep);
%save data to object structure
linsys.taylor.eAt=eAt;

F=linsys.taylor.F;
RV=linsys.taylor.RV;
inputCorr=linsys.taylor.inputCorr;
if iscell(linsys.taylor.Rtrans)
    Rtrans=linsys.taylor.Rtrans{1};
else
    Rtrans=linsys.taylor.Rtrans;
end


%first time step homogeneous solution
Rhom_tp=eAt*Rinit + Rtrans;
if isa(Rinit,'polyZonotope') || isa(Rinit,'conPolyZono')
    Rhom=enclose(Rinit,Rhom_tp)+F*zonotope(Rinit)+inputCorr;
elseif isa(Rinit,'zonoBundle') 
    Rhom=enclose(Rinit,Rhom_tp)+F*Rinit.Z{1}+inputCorr;
else
    try
        Rhom=enclose(Rinit,Rhom_tp)+F*Rinit+inputCorr;
    catch
        Rhom=enclose(Rinit,Rhom_tp)+F*interval(Rinit)+inputCorr;
    end
end

%reduce zonotopes
Rhom = reduce(Rhom,options.reductionTechnique,options.zonotopeOrder);
Rhom_tp = reduce(Rhom_tp,options.reductionTechnique,options.zonotopeOrder);
if ~isnumeric(RV)
    RV = reduce(RV,options.reductionTechnique,options.zonotopeOrder);
end

%save homogeneous and particulate solution
options.Rhom=Rhom;
options.Rhom_tp=Rhom_tp;
options.Raux=RV;
if strcmp(options.linAlg,'wrapping-free')
    options.Rpar=interval(RV);
else
    options.Rpar=RV;
end
options.Rtrans=linsys.taylor.Rtrans;

%total solution
if isa(Rinit,'polytope')
    %convert zonotopes to polytopes
    Radd=polytope(RV);
    Rtotal=Rhom+Radd;
    Rtotal_tp=Rhom_tp+Radd;
else
    %original computation
    Rtotal=Rhom+RV;
    Rtotal_tp=Rhom_tp+RV;
end

%write results to reachable set struct Rfirst
Rfirst.tp=Rtotal_tp;
Rfirst.ti=Rtotal;

end


% Auxiliary functions -----------------------------------------------------

function obj = aux_exponential(obj,options)
% computes the overapproximation of the exponential of a system matrix up
% up to a certain accuracy
   
%load data from object/options structure
A=obj.A;
A_abs=abs(A);
taylorTerms=options.taylorTerms;
n=obj.nrOfDims;
factors = options.factor;

%initialize
Apower = cell(taylorTerms,1);
Apower{1} = A;  
Apower_abs = cell(taylorTerms,1);
Apower_abs{1} = A_abs; 
M = eye(n);
    
%compute powers for each term and sum of these
for i=1:taylorTerms
    %compute powers
    Apower{i+1}=Apower{i}*A;
    Apower_abs{i+1}=Apower_abs{i}*A_abs;
    M = M + Apower_abs{i}*factors(i);
end   
%determine error due to finite Taylor series, see Prop. (2) in [1]
W=expm(A_abs*options.timeStep)-M;
%compute absolute value of W for numerical stability
W=abs(W);
E=interval(-W,W);
    
%write to object structure
obj.taylor.powers=Apower;
obj.taylor.error=E;

end

function obj = aux_tie(obj,options)
% computes the error done by building the convex hull of time point solutions

%load data from object/options structure
Apower=obj.taylor.powers;
taylorTerms=options.taylorTerms;
rbyfac=options.factor;
n=obj.nrOfDims;

%initialize Asum
Asum_pos=zeros(n);
Asum_neg=zeros(n);

for i=2:taylorTerms
    %compute factor
    exp1=-i/(i-1); exp2=-1/(i-1);
    factor=(i^exp1-i^exp2)*rbyfac(i); 
    
    %init Apos, Aneg
    Apos=zeros(n);
    Aneg=zeros(n);
    
    %obtain positive and negative parts
    pos_ind = Apower{i}>0;
    neg_ind = Apower{i}<0;
    
    Apos(pos_ind) = Apower{i}(pos_ind);
    Aneg(neg_ind) = Apower{i}(neg_ind);
    
    %compute powers; factor is always negative
    Asum_pos=Asum_pos + factor*Aneg; 
    Asum_neg=Asum_neg + factor*Apos;
end
%instantiate interval matrix
Asum = interval(Asum_neg,Asum_pos);

%write to object structure
obj.taylor.F=Asum+obj.taylor.error;

end

function [obj,options] = aux_inputSolution(obj,params,options)
% computes the bloating due to the input

% possible effect of disturbance: total input is then
%   V = obj.B*options.U + (obj.E*options.W - center(obj.E*options.W))
%   vTrans = obj.B*options.uTrans + center(obj.E*options.W)
if isfield(params,'W')
    params.W = obj.E*params.W;
    Wcenter = center(params.W);
    W = params.W + (-Wcenter);
else % entered by linearized system, e.g., from nonlinearSys
    Wcenter = 0; W = 0;
end

% set of possible inputs
V = obj.B*params.U + W;
% compute vTrans (including disturbance center and constant offset)
vTrans = obj.B*params.uTrans + Wcenter + obj.c;

% do we have a time-varying input solution?
options.isRV = ~representsa_(V,'origin',eps);

Apower = obj.taylor.powers;
E = obj.taylor.error;
r = options.timeStep;
n = obj.nrOfDims;
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
    obj = aux_inputTie(obj,options);
    inputF = obj.taylor.inputF;
    inputCorr = inputF*zonotope(vTrans);
end


%write to object structure
obj.taylor.V = V;
obj.taylor.RV = inputSolV;
obj.taylor.Rtrans = inputSolVtrans;
obj.taylor.inputCorr = inputCorr;
obj.taylor.eAtInt = eAtInt;

end

function obj = aux_inputTie(obj,options)
% computes the error done by the linear assumption of the constant input
% solution

% load data from object structure
Apower=obj.taylor.powers;
E = obj.taylor.error;
taylorTerms=options.taylorTerms;
r=options.timeStep;
n=obj.nrOfDims;

% initialize Asum
Asum_pos=zeros(n);
Asum_neg=zeros(n);

for i=2:(taylorTerms+1)
    % compute factor
    exp1=-i/(i-1); exp2=-1/(i-1);
    factor=(i^exp1-i^exp2)*options.factor(i); 
    
    % init Apos, Aneg
    Apos=zeros(n);
    Aneg=zeros(n);
    
    % obtain positive and negative parts
    pos_ind = Apower{i-1}>0;
    neg_ind = Apower{i-1}<0;
    
    Apos(pos_ind) = Apower{i-1}(pos_ind);
    Aneg(neg_ind) = Apower{i-1}(neg_ind);
    
    % compute powers; factor is always negative
    Asum_pos=Asum_pos + factor*Aneg; 
    Asum_neg=Asum_neg + factor*Apos;
end
% instantiate interval matrix
Asum = interval(Asum_neg,Asum_pos);

% compute error due to finite Taylor series according to internal document
% "Input Error Bounds in Reachability Analysis"
Einput = E*r;

%write to object structure
obj.taylor.inputF=Asum+Einput; %rewrite this equation when E is computed with the new method

end

% ------------------------------ END OF CODE ------------------------------
