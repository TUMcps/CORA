function [Rti,Rtp,dimForSplit,options] = linReach(sys,Rstart,params,options)
% linReach - computes the reachable set after linearization
%
% Syntax:
%    [Rti,Rtp,dimForSplit,options] = linReach(sys,Rstart,params,options)
%
% Inputs:
%    sys - nonlinearSys or nonlinParamSys object
%    Rstart - initial reachable set
%    params - model parameters
%    options - struct with algorithm settings
%
% Outputs:
%    Rti - reachable set for time interval
%    Rtp - reachable set for time point
%    dimForSplit - dimension that is split to reduce the lin. error
%    options - struct with algorithm settings
%
% References: 
%   [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%       uncertain parameters using conservative linearization"
%   [2] M. Althoff. "Reachability analysis of nonlinear systems using 
%       conservative polynomialization and non-convex sets"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: initReach, post

% Authors:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       17-January-2008
% Last update:   29-June-2009
%                23-July-2009
%                10-July-2012
%                18-September-2012
%                09-August-2016
%                12-September-2017
%                02-January-2020 (NK, restructured the function)
%                22-April-2020 (MW, simplification)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% extract initial set and abstraction error
Rinit = Rstart.set;
abstrerr = Rstart.error;

% necessary to update part of abstraction that is dependent on x0 when
% linearization remainder is not computed
if isfield(options,'updateInitFnc')
    currentStep = round((options.t-params.tStart)/options.timeStep)+1;
    Rinit = options.updateInitFnc(Rinit,currentStep);
end

% linearize the nonlinear system
[sys,linsys,linParams,linOptions] = linearize(sys,Rinit,params,options); 

% translate Rinit by linearization point
Rdelta = Rinit + (-sys.linError.p.x);

% compute reachable set of the linearized system
if isa(sys,'nonlinParamSys') && isa(params.paramInt,'interval')
    [linsys,R] = initReach_inputDependence(linsys,Rdelta,linParams,linOptions);
    Rtp = R.tp; Rti = R.ti;
elseif isa(linsys,'linParamSys')
    R = initReach(linsys,Rdelta,linParams,linOptions);
    Rtp = R.tp; Rti = R.ti;
else
    [Rtp,Rti,~,~,PU,Pu,~,C_input] = oneStep(linsys,Rdelta,...
        linParams.U,linParams.uTrans,options.timeStep,options.taylorTerms);

    if strcmp(options.alg,'poly')
        % pre-compute set of state differences
        Rdiff = aux_deltaReach(linsys,Rdelta,PU,Pu,C_input,...
            options.timeStep,options.taylorTerms,...
            options.reductionTechnique,options.intermediateOrder);
        
        % pre-compute static abstraction error
        if options.tensorOrder > 2
            [H,Zdelta,errorStat,T,ind3,Zdelta3] = ...
                priv_precompStatError(sys,Rdelta,params,options);
        end
    end
end
if isfield(options,'approxDepOnly') && options.approxDepOnly
    if ~exist('errorStat','var')
        errorStat = [];
    end
    R.tp = Rtp; R.ti = Rti;
    [Rtp,Rti,dimForSplit,options] = aux_approxDepReachOnly(linsys,sys,R,options,errorStat);
    return;
end
% compute reachable set of the abstracted system including the
% abstraction error using the selected algorithm
if strcmp(options.alg,'linRem')
    [Rtp,Rti,perfInd] = aux_linReach_linRem(sys,R,Rinit,Rdelta,params,options);
else

    % loop until the actual abstraction error is smaller than the 
    % estimated linearization error
    perfIndCurr = Inf; perfInd = 0;
    
    % used in AROC for reachsetOptimalControl (reachSet with previous
    % linearization error)
    if isfield(options,'prevErr')
        perfIndCurr = 1;
        if isfield(options,'prevErrScale')
            scale = options.prevErrScale;
        else 
            scale = 0.8;
        end
        Rerror = scale*options.prevErr;
    end
    
    while perfIndCurr > 1 && perfInd <= 1
        % estimate the abstraction error 
        appliedError = 1.1*abstrerr;
        Verror = zonotope(0*appliedError,diag(appliedError));
        if isa(linsys,'linParamSys')
            RallError = errorSolution(linsys,options,Verror);
        else
            RallError = particularSolution_timeVarying(linsys,...
                Verror,options.timeStep,options.taylorTerms);
        end

        % compute the abstraction error using the conservative
        % linearization approach described in [1]
        if strcmp(options.alg,'lin')

            % compute overall reachable set including linearization error
            Rmax = Rti+RallError;
            % compute linearization error
            [trueError,VerrorDyn] = priv_abstrerr_lin(sys,Rmax,params,options);
            VerrorStat = zeros(sys.nrOfDims,1);
            
        % compute the abstraction error using the conservative
        % polynomialization approach described in [2]    
        else

            % compute overall reachable set including linearization error
            Rmax = Rdelta+zonotope(Rdiff)+RallError;
            % compute abstraction error
            [trueError,VerrorDyn,VerrorStat] = ...
                priv_abstrerr_poly(sys,Rmax,Rdiff+RallError,params,options, ...
                                    H,Zdelta,errorStat,T,ind3,Zdelta3);

        end

        % compare linearization error with the maximum allowed error
        perfIndCurr = max(trueError./appliedError);    
        perfInd = max(trueError./options.maxError);

        abstrerr = trueError;
        
        % clean exit in case of set explosion
        if any(abstrerr > 1e+100)
            throw(CORAerror('CORA:reachSetExplosion'));
        end
    end

    % translate reachable sets by linearization point
    Rti = Rti+sys.linError.p.x;
    Rtp = Rtp+sys.linError.p.x;

    % compute the reachable set due to the linearization error
    if ~exist('Rerror','var')
        if isa(linsys,'linParamSys')
            Rerror = errorSolution(linsys,options,VerrorDyn);
        else
            Rerror_dyn = particularSolution_timeVarying(linsys,...
                VerrorDyn,options.timeStep,options.taylorTerms);
            Rerror_stat = particularSolution_constant(linsys,...
                VerrorStat,options.timeStep,options.taylorTerms);
            Rerror = Rerror_dyn + Rerror_stat;
        end
        if isfield(options,'approxErr') && options.approxErr
            options.prevErr = Rerror;
        end
    end
    
    % add the abstraction error to the reachable sets
    if strcmp(options.alg,'poly') && (isa(Rerror,'polyZonotope') || ...
                                      isa(Rerror,'conPolyZono'))
        Rti=exactPlus(Rti,Rerror);
        Rtp=exactPlus(Rtp,Rerror);
    else
        Rti=Rti+Rerror;
        Rtp=Rtp+Rerror;
    end
end

% determine the best dimension to split the set in order to reduce the
% linearization error
dimForSplit = [];

if perfInd > 1
    dimForSplit = priv_select(sys,Rstart,params,options);
end

% store the linearization error
Rtp_.set = Rtp;
Rtp_.error = abstrerr;
Rtp = Rtp_;


end


% Auxiliary functions -----------------------------------------------------

function Rdelta = aux_deltaReach(sys,Rinit,RV,Rtrans,inputCorr,...
    timeStep,truncationOrder,reductionTechnique,intermediateOrder)
% compute delta reachable set (set of states differences): the notable
% different to linearSys functions is that we have (e^At - I)*Rinit and an
% enclose-call with the origin

% compute/read out helper variables
options = struct('timeStep',timeStep,'ithpower',truncationOrder);
eAt = getTaylor(sys,'eAdt',options);
F = readFieldForTimeStep(sys.taylor,'F',timeStep);

% first time step homogeneous solution
n = sys.nrOfDims;
Rhom_tp_delta = (eAt - eye(n))*Rinit + Rtrans;

if isa(Rinit,'zonotope')
    %original computation
    O = zonotope.origin(n);
    Rhom=enclose(O,Rhom_tp_delta)+F*Rinit+inputCorr;
elseif isa(Rinit,'polyZonotope') || isa(Rinit,'conPolyZono')
    O = zeros(n)*Rhom_tp_delta;  % to retain dependencies!
    Rhom=enclose(O,Rhom_tp_delta)+F*zonotope(Rinit)+inputCorr;
elseif isa(Rinit,'zonoBundle')
    O = zonoBundle.origin(n);
    Rhom=enclose(O,Rhom_tp_delta)+F*Rinit.Z{1}+inputCorr;
end

%reduce zonotope
Rhom=reduce(Rhom,reductionTechnique,intermediateOrder);
if ~isnumeric(RV)
    RV=reduce(RV,reductionTechnique,intermediateOrder);
end

%total solution
if isa(Rinit,'polytope')
    %convert zonotopes to polytopes
    Radd=polytope(RV);
    Rdelta=Rhom+Radd;
else
    %original computation
    Rdelta=Rhom+RV;
end

end

function [Rtp,Rti,perfInd] = aux_linReach_linRem(sys,R,Rinit,Rdelta,params,options)
% Compute the reachable set for the linearized system using an algorithm
% that is based on the linearization of the Lagrange remainder
    
    % compute the reachable set for the linearized system
    options.alg = 'lin';
    
    [sys,linsys,linParams,linOptions] = linearize(sys,Rinit,params,options);
    if isa(linsys,'linParamSys')
        linOptions.compTimePoint = true;
    end
    if isa(sys,'nonlinParamSys') && isa(params.paramInt,'interval')
        [~,Rlin] = initReach_inputDependence(linsys,Rdelta,linParams,linOptions);
        Ro_int = interval(Rlin.ti);
    elseif isa(linsys,'linParamSys')
        Rlin = initReach(linsys,Rdelta,linParams,linOptions);
        Ro_int = interval(Rlin.ti);
    elseif isa(linsys,'linearSys')
        Rlinti = oneStep(linsys,Rdelta,linParams.U,...
            linParams.uTrans,options.timeStep,options.taylorTerms);
        Ro_int = interval(Rlinti);
    end
    
    % compare the computed reachable set to the reachable set of the
    % linearized system in order to decide if splitting is required
    Rti_int = interval(R.ti);
    
    trueError = max(abs(Rti_int.inf-Ro_int.inf),abs(Rti_int.sup-Ro_int.sup));
    perfInd = max(trueError./options.maxError);
    
    % translate reachable sets by linearization point
    Rti = R.ti + sys.linError.p.x;
    Rtp = R.tp + sys.linError.p.x;
    
end

% TODO: put this somewhere else
function [Rtp,Rti,dimForSplit,options] = aux_approxDepReachOnly(linsys,nlnsys,R,options,errorStat)
% Computes an approximation of the reachable set for controller synthesis.
% Compared to the over-approximative reachability algorithm, the
% higher-order terms are only evaluated for the time-point reachable set
% (errorStat) and the Lagrange remainder is neglected.

    %read tp and ti
    R_tp = R.tp;
    R_ti = R.ti;

    if representsa_(errorStat,'emptySet',eps) || all(representsa_(errorStat,'origin',eps))
        % we do not need to consider errorStat then
        R_tp = R_tp + nlnsys.linError.p.x;
        R_ti = R_ti + nlnsys.linError.p.x;
    else
        % consider errorStat
        [id,~,ind] = unique(R_ti.id); E = zeros(length(id),size(R_ti.E,2));
        for i = 1:length(ind)
            E(ind(i),:) = E(ind(i),:) + R_ti.E(i,:);
        end
        R_ti = polyZonotope(R_ti.c,R_ti.G,R_ti.GI,E,id);

        Asum = options.timeStep*eye(linsys.nrOfDims);
        for i = 1:options.taylorTerms
            Asum = Asum + linsys.taylor.Apower{i}*linsys.taylor.dtoverfac{1}(i+1);
        end
        eAtInt = Asum + linsys.taylor.E{1}*options.timeStep;

        Rerror = eAtInt*errorStat;
        R_tp = exactPlus(R_tp,Rerror) + nlnsys.linError.p.x;
        R_ti = exactPlus(R_ti,Rerror) + nlnsys.linError.p.x;
    end

   
    % init output variables
    R_tp = noIndep(reduce(R_tp,options.reductionTechnique,options.zonotopeOrder));
    R_ti = noIndep(reduce(R_ti,options.reductionTechnique,options.zonotopeOrder));
    Rtp_.set = R_tp;
    Rti = R_ti;
    Rtp_.error = zeros(length(R_tp.c),1);
    Rtp = Rtp_;
    dimForSplit = [];
end

% ------------------------------ END OF CODE ------------------------------
