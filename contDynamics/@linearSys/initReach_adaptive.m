function [Rend,options,linOptions] = initReach_adaptive(linsys,Rstart,options,linParams,linOptions)
% initReach_adaptive - computes the linearized reachable set
%    from an originally nonlinear system (linearization error = 0),
%    loops until time step found which satisfies abstraction error bound
%
% Syntax:
%    Rend = initReach_adaptive(linsys,Rstart,options,linParams,linOptions)
%
% Inputs:
%    linsys - linearized system
%    Rstart - reachable set of current time point
%    options - options for nonlinear system
%    linParams - model parameter for linearized system
%    linOptions - options for linearized system
%
% Outputs:
%    Rend - reachable set (time interval/point) of current time + time step
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linReach

% Authors:       Mark Wetzlinger
% Written:       25-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% re-use same truncation order as in last iteration of same step to
% facilitate reduction by indices (particular solutions must have same
% number of generators across iterations in the same step)
if isfield(options,'tt_lin') && length(options.tt_lin) == options.i
    linOptions.taylorTerms = options.tt_lin(options.i);
end
% exponential matrix and time interval error (incl. adaptive taylorTerms)
[linsys,linOptions] = aux_expmtie_adaptive(linsys,linOptions);
% compute reachable set due to input
[linsys,linOptions] = aux_inputSolution(linsys,linParams,linOptions);
%change the time step
linsys.taylor.timeStep = linOptions.timeStep;
% save taylorTerms for next iteration of same step
options.tt_lin(options.i,1) = linOptions.taylorTerms;
options.etalinFro(options.i,1) = linOptions.etalinFro;

%compute reachable set of first time interval
eAt = expm(linsys.A*linOptions.timeStep);
%save data to object structure
linsys.taylor.eAt=eAt;

F=linsys.taylor.F;
inputCorr=linsys.taylor.inputCorr;
Rtrans=linsys.taylor.Rtrans;

%first time step homogeneous solution
Rhom_tp = eAt*Rstart + Rtrans;
Rhom = enclose(Rstart,Rhom_tp) + F*Rstart + inputCorr;

% preliminary solutions without RV
if isfield(options,'gredIdx')
    if length(options.gredIdx.Rhomti) == options.i
        % select the same generators for reduction as in the previous
        % iteration of the current time step to limit the influence of
        % the order reduction on the computation of the gain order
        Rend.ti = reduce(Rhom,'idx',options.gredIdx.Rhomti{options.i});
    else
        % reduce adaptively and store indices of generators that have
        % not been selected for reduction
        [Rend.ti,~,options.gredIdx.Rhomti{options.i}] = ...
            reduce(Rhom,'adaptive',options.redFactor);
    end
else
    % safety branch (call from algorithm without gredIdx)
    Rend.ti = reduce(Rhom,'adaptive',options.redFactor);
end

if isfield(options,'gredIdx')
    if length(options.gredIdx.Rhomtp) == options.i
        % select the same generators for reduction as in the previous
        % iteration of the current time step to limit the influence of
        % the order reduction on the computation of the gain order
        Rend.tp = reduce(Rhom_tp,'idx',options.gredIdx.Rhomtp{options.i});
    else
        % reduce adaptively and store indices of generators that have
        % not been selected for reduction
        [Rend.tp,~,options.gredIdx.Rhomtp{options.i}] = ...
            reduce(Rhom_tp,'adaptive',options.redFactor);
    end
else
    % safety branch (call from algorithm without gredIdx)
    Rend.tp = reduce(Rhom_tp,'adaptive',options.redFactor);
end

% reduce and add RV only if exists
if linOptions.isRV
    % read and reduce RV from struct
    if isfield(options,'gredIdx')
        if length(options.gredIdx.Rpar) == options.i
            % select the same generators for reduction as in the previous
            % iteration of the current time step to limit the influence of
            % the order reduction on the computation of the gain order
            RV = reduce(linsys.taylor.RV,'idx',options.gredIdx.Rpar{options.i});
        else
            % reduce adaptively and store indices of generators that have
            % not been selected for reduction
            [RV,~,options.gredIdx.Rpar{options.i}] = ...
                reduce(linsys.taylor.RV,'adaptive',options.redFactor);
        end
    else
        % safety branch (call from algorithm without gredIdx)
        RV = reduce(linsys.taylor.RV,'adaptive',options.redFactor);
    end

    %total solution
    Rend.ti = Rend.ti + RV;
    Rend.tp = Rend.tp + RV;
end

end


% Auxiliary functions -----------------------------------------------------

function [linsys,options] = aux_expmtie_adaptive(linsys,options)
% computes the remainder of the exponential matrix and the correction matrix

% load data from object/options structure
A = linsys.A;
A_abs = abs(A);
n = linsys.nrOfDims;
deltat = options.timeStep;
% deltatbyfac = options.factor;
taylorTermsGiven = isfield(options,'taylorTerms');

% initialize 
Apower{1} = A;
Apower_abs{1} = A_abs;
M = eye(n);
Asum_pos = zeros(n);
Asum_neg = zeros(n);
    
% increment taylorTerms until convergence
eta = 1;
while true
    % exponential: 1:eta
    % compute powers
    Apower{eta+1} = Apower{eta}*A;
    Apower_abs{eta+1} = Apower_abs{eta}*A_abs;
    deltatbyfac(eta,1) = deltat^eta / factorial(eta);
    M = M + Apower_abs{eta}*deltatbyfac(eta);
    
    if eta==1; eta = eta+1; continue; end % F starts at eta=2
    
    % tie: 2:eta
    % compute factor
    exp1 = -(eta)/(eta-1); exp2 = -1/(eta-1);
    factor = ((eta)^exp1-(eta)^exp2) * deltatbyfac(eta); 
    
    % init Apos, Aneg
    Apos = zeros(n);
    Aneg = zeros(n);
    
    % obtain positive and negative parts
    pos_ind = Apower{eta} > 0;
    neg_ind = Apower{eta} < 0;
    
    Apos(pos_ind) = Apower{eta}(pos_ind);
    Aneg(neg_ind) = Apower{eta}(neg_ind);
    
    % compute powers; factor is always negative
    Asum_pos = Asum_pos + factor*Aneg; 
    Asum_neg = Asum_neg + factor*Apos;
    
    % norm as stopping criterion (no abs() because already all positive)
    normAfro(eta,1) = norm(Asum_pos - Asum_neg,'fro');

    if taylorTermsGiven
        % has to be same eta as in other runs of same step
        stopCondition = eta == options.taylorTerms;
    else
        stopCondition = ~any(any(Apower{eta})) ... % nilpotent
            || 1 - normAfro(eta-1)/normAfro(eta) < options.zetaTlin; % relative convergence
    end

    if stopCondition
        % determine error due to finite Taylor series, see eq.(2) in [1]
        W = expm(A_abs*options.timeStep) - M;
        % compute absolute value of W for numerical stability
        W = abs(W);
        E = interval(-W,W);
        % instantiate interval matrix
        Asum = interval(Asum_neg,Asum_pos);
        deltatbyfac(eta+1,1) = deltat^(eta+1) / factorial(eta+1); %inputSolution
        options.taylorTerms = eta;
        options.factor = deltatbyfac;
        options.etalinFro = 1 - normAfro(eta-1)/normAfro(eta);
        break
    else
        eta = eta + 1;
        if eta > 50
            throw(MException('expmtie:notconverging',...
                'Time Step Size too big.'));
        end
    end
end

% write to object structure
linsys.taylor.powers = Apower;
linsys.taylor.error = E;
linsys.taylor.F = Asum + E;

end

function [linsys,options] = aux_inputSolution(linsys,params,options)
% computes the bloating due to the input

% possible effect of disturbance: total input is then
%   V = obj.B*options.U + (obj.E*options.W - center(obj.E*options.W))
%   vTrans = obj.B*options.uTrans + center(obj.E*options.W)
if isfield(params,'W')
    params.W = linsys.E*params.W;
    Wcenter = center(params.W);
    W = params.W + (-Wcenter);
else % entered by linearized system, e.g., from nonlinearSys
    Wcenter = 0; W = 0;
end

% set of possible inputs
V = linsys.B*params.U + W;
% compute vTrans (including disturbance center and constant offset)
vTrans = linsys.B*params.uTrans + Wcenter + linsys.c;

% do we have a time-varying input solution?
options.isRV = ~representsa_(V,'origin',eps);

Apower = linsys.taylor.powers;
E = linsys.taylor.error;
r = options.timeStep;
n = linsys.nrOfDims;
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
    linsys = aux_inputTie(linsys,options);
    inputF = linsys.taylor.inputF;
    inputCorr = inputF*zonotope(vTrans);
end


%write to object structure
linsys.taylor.V = V;
linsys.taylor.RV = inputSolV;
linsys.taylor.Rtrans = inputSolVtrans;
linsys.taylor.inputCorr = inputCorr;
linsys.taylor.eAtInt = eAtInt;

end

function linsys = aux_inputTie(linsys,options)
% computes the error done by the linear assumption of the constant input
% solution

% load data from object structure
Apower=linsys.taylor.powers;
E = linsys.taylor.error;
taylorTerms=options.taylorTerms;
r=options.timeStep;
n=linsys.nrOfDims;

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
linsys.taylor.inputF=Asum+Einput; %rewrite this equation when E is computed with the new method

end

% ------------------------------ END OF CODE ------------------------------
