function [params,options] = postProcessing(sys,func,params,options)
% postProcessing - perform post-processing operations for
%    model parameters and algorithm parameters
%
% Syntax:
%    [params,options] = postProcessing(sys,func,params,options)
%
% Inputs:
%    sys - contDynamics or hybridDynamics object
%    func - function name ('reach', 'simulateRandom', etc.)
%    params - user-defined model parameters
%    options - user-defined algorithm parameters
%
% Outputs:
%    params - user-defined model parameters updated for internal usage
%    options - user-defined algorithm parameters updated for internal usage
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      04-February-2021
% Last update:  07-July-2021 (MA, U should not be converted for
%                                 simulateGaussian and observe)
% Last revision:---

%------------- BEGIN CODE --------------


% only contDynamics -------------------------------------------------------

if isa(sys,'contDynamics')
    
    % convert U to a zonotope if given as an interval
    if ~any(strcmp(func,{'reachInnerProjection','simulateRandom','observe'}))
        [params,options] = convert_U(sys,params,options);
    end
    
    % set linAlg option for nonlinear system classes (calling linSys)
    if ~isa(sys,'linearSys') && ~isa(sys,'linearSysDT')
        [params,options] = set_linAlg(sys,params,options);
    end
    
    % extend vectors to correct length if default setting
    if strcmp(func,'observe')
        [params,options] = set_u_y(sys,params,options);
    end
    
    % set input set U and input vector uTrans|uTransVec
    if ~(contains(func,'reachInner') || strcmp(func,'observe'))
        
        % correct splitting of inputs into U (containing origin) and uTrans
        % (if shift is constant) or uTransVec (if shift is time-varying)
        [params,options] = set_U_uTrans_uTransVec(sys,params,options);
        % adjust params.uTransVec|tu for verify
        if contains(func,'verify')
            [params,options] = set_u_tu(sys,params,options);
        end
    
        if strcmp(func,'simulateRandom')
            % set u, tu, and nrConstInp for simulation
            [params,options] = set_u_nrConstInp_tu(sys,params,options);
        elseif (isa(sys,'linearSys') || isa(sys,'linParamSys') || isa(sys,'linProbSys'))
            % set internal value originContained
            [params,options] = set_originContained(sys,params,options);
        end
        
    end
    
    % set alg and N
    if isa(sys,'nonlinearSys') && contains(func,'reachInner')
        if strcmp(func,'reachInnerScaling')
            % use polynomialization, compute number of steps
            [params,options] = set_R0_alg_tensorOrder(sys,params,options);
            [params,options] = set_N(sys,params,options);
        elseif strcmp(func,'reachInnerProjection')
            % R0 converted to an interval
            [params,options] = set_R0int(sys,params,options);
        end
    end
    
    % set polyZonotope restructuring
    if ~strcmp(func,'simulateRandom')
        if isa(sys,'nonlinearSysDT')
            [params,options] = set_volApproxMethod(sys,params,options);
        elseif isa(sys,'nonlinearSys') || isa(sys,'nonlinParamSys')
            [params,options] = set_volApproxMethod(sys,params,options);
        end
    end
    
    % nonlinear adaptive tuning
    if strcmp(func,'reach') && contains(class(sys),'nonlin') && ...
            contains(options.alg,'adaptive')
        if isa(sys,'nonlinearSys')
            [params,options] = set_nonlinearSys_adaptive(sys,params,options);
        elseif isa(sys,'nonlinearSysDT')
            [params,options] = set_nonlinearSysDT_adaptive(sys,params,options);
        elseif isa(sys,'nonlinDASys')
            [params,options] = set_nonlinDASys_adaptive(sys,params,options);
        end
    end

    % verification: specifications and corresponding time intervals
    if strcmp(func,'verify')
        % initialize time intervals where (un)safe sets are not yet verified
        [params,options] = set_unsatIntervals(sys,params,options);

        % check if computation of inner-approximation can be skipped in favor
	    % of a simpler method (only if F/G consists of only of halfspace)
        [params,options] = set_safeUnsafeSetFastInner(sys,params,options);
    end
    
end

% only hybridAutomaton ----------------------------------------------------

if isa(sys,'hybridAutomaton')
   
    if strcmp(func,'reach')
        % set timeStep for each location
        [params,options] = set_HA_timeStep(sys,params,options);
    end
    if any(strcmp(func,{'reach','simulateRandom'}))
        % set input for locations
        [params,options] = set_HA_inputs(sys,params,options);
    end
    if strcmp(func,'simulate')
        % set input for locations
        [params,options] = set_HA_sim_inputs(sys,params,options);
    end
    
end

% only parallelHybridAutomaton --------------------------------------------

if isa(sys,'parallelHybridAutomaton')
   
    % set input for locations
    if any(strcmp(func,{'reach','simulateRandom'}))
        [params,options] = set_pHA_inputs(sys,params,options);
    end
    if strcmp(func,'simulate')
        [params,options] = set_pHA_sim_inputs(sys,params,options);
    end
    
end

end



% Auxiliary Functions: contDynamics

function [params,options] = convert_U(sys,params,options)

% convert U to a zonotope if given as an interval
if isa(params.U,'interval')
    params.U = zonotope(params.U);
end

end

function [params,options] = set_linAlg(sys,params,options)

% only 'standard' for all linearized systems
options.linAlg = 'standard';

end

function [params,options] = set_U_uTrans_uTransVec(sys,params,options)

% shift u by center of U
centerU = center(params.U);
if any(centerU)
    params.u = params.u + center(params.U);
    params.U = params.U + (-centerU);
end

% write u to uTrans or uTransVec
inputTrajLength = size(params.u,2);
if inputTrajLength > 1
    params.uTransVec = params.u;
else
    params.uTrans = params.u;
end
params = rmfield(params,'u');

end

function [params,options] = set_u_tu(sys,params,options)
% required: params.uTransVec, params.tu
% removes duplicates and adjusts params.uTransVec

% collapse uTransVec if necessary
if isfield(params,'uTransVec') && isfield(params,'tu')
    idxChange = [1;1+find(any(diff(params.uTransVec'),2))];
    if length(idxChange) <= size(params.uTransVec,2)
        params.uTransVec = params.uTransVec(:,idxChange);
        params.tu = params.tu(idxChange);
    end
end

end

function [params,options] = set_u_nrConstInp_tu(sys,params,options)

if ~strcmp(options.type,'rrt')
    
    % adapt options.nrConstInp: make vector
    if isscalar(options.nrConstInp)
        if ~isfield(params,'uTransVec') && isscalar(params.tu) ...
                && withinTol(params.tu(1),params.tStart)
            stepsize = (params.tFinal - params.tStart) / options.nrConstInp;
            params.tu = params.tStart:stepsize:params.tFinal;
        end
        options.nrConstInp = repelem(1,options.nrConstInp);
    end

    % append final time to params.tu unless already provided
    if abs(params.tu(end) - params.tFinal) > 1e-12
        try
            params.tu = [params.tu; params.tFinal];
        catch
            params.tu = [params.tu, params.tFinal];
        end
    end

end

end

function [params,options] = set_originContained(sys,params,options)

U = params.U;
% define actual input set
if isfield(params,'uTrans')
    U = U + params.uTrans;
elseif isfield(params,'uTransVec')
    % input trajectory
    options.originContained = false; return;
end

% transform scalar identity matrix to matrix
B = sys.B;
if isscalar(B)
    if B == 1
        B = eye(sys.dim);
    elseif B == 0
        B = zeros(sys.dim);
    end
end

% options.originContained = false;
if isa(sys,'linearSys') && ~isempty(sys.c)
    % check if any of equality constraints B(i)*U = -c(i) = 0
    % is unsatisfiable
    found = false;
    for i = 1:size(B,1)
       hp = conHyperplane(B(i,:),-sys.c(i));
       if ~isIntersecting(hp,U)
          options.originContained = false;
          found = true;
       end
    end
    if ~found
        % check if vTrans = B*U + c contains the origin
        vTrans = B*U + sys.c;
        % faster computation if vTrans is an interval
        if isInterval(vTrans)
            vTransInt = interval(vTrans);
            options.originContained = contains(vTransInt,zeros(dim(vTrans),1));
        else
            if isa(vTrans,'polyZonotope')
                options.originContained = contains(vTrans,zeros(dim(vTrans),1),'approx');
            else
                options.originContained = contains(vTrans,zeros(dim(vTrans),1));
            end
        end
    end
elseif isInterval(U)
    % no constant input, faster computation if U is an interval
    I = interval(U);
    options.originContained = contains(I,zeros(dim(U),1));
else
    % no constant input c, U not an interval
    options.originContained = contains(U,zeros(dim(U),1));
end

end

function [params,options] = set_u_y(sys,params,options)

if size(params.u,2) == 1 && ~any(params.u)
    % extend vector to correct number of columns
    params.u = repmat(params.u,1,(params.tFinal-params.tStart)/sys.dt);   
end
% rewrite to uTransVec
params.uTransVec = params.u;
params = rmfield(params,'u');

if size(params.y,2) == 1 && ~any(params.y)
    % extend vector to correct number of columns
    params.y = repmat(params.y,1,(params.tFinal-params.tStart)/sys.dt);   
end

end

function [params,options] = set_R0_alg_tensorOrder(sys,params,options)

params.R0 = polyZonotope(params.R0);
options.alg = 'poly';
options.tensorOrder = 3;

end

function [params,options] = set_R0int(sys,params,options)

params.R0 = interval(params.R0);

end

function [params,options] = set_volApproxMethod(sys,params,options)

if isa(sys,'nonlinearSysDT')
    if isFullDim(params.R0)
        options.polyZono.volApproxMethod = 'interval';
    else
        options.polyZono.volApproxMethod = 'pca';
    end
else
    % isfield because of nonlinearSys/reachInner
    if isfield(options,'alg') && strcmp(options.alg,'poly')
        if isFullDim(params.R0)
            options.polyZono.volApproxMethod = 'interval';
        else
            options.polyZono.volApproxMethod = 'pca';
        end
    end
end

end

function [params,options] = set_N(sys,params,options)

if isfield(options,'timeStepInner')
    options.N = round(options.timeStepInner/options.timeStep,0);
    options = rmfield(options,'timeStepInner');
else
    options.N = length(params.tStart:options.timeStep:params.tFinal)-1;
end

end

function [params,options] = set_nonlinearSys_adaptive(sys,params,options)

options.decrFactor = 0.90;              % adaptation of delta t
options.zetaTlin = 0.0005;              % zeta_T,lin (taylorTerms)
options.zetaTabs = 0.005;               % zeta_T,abs (taylorTerms)
if contains(options.alg,'lin')
    params.R0 = zonotope(params.R0);    % convert to zono
    params.R = params.R0;               % for start set
    options.redFactor = 0.0005;         % zeta_Z (zonotope order)
    options.zetaK = 0.90;               % zeta_K (tensorOrder)
    % zeta_h (timeStep) ... depends on order (chosen in code) and decrFactor
    options.zetaphi = [0.85; 0.76; 0.68]; 	
elseif contains(options.alg,'poly')
    % options for poly (fixed in adaptive)
    options.polyZono.volApproxMethod = 'interval';           % polyZono
    options.polyZono.maxDepGenOrder = 50;                    % polyZono (unused...)
    options.polyZono.maxPolyZonoRatio = 0.05;                % polyZono
    options.polyZono.restructureTechnique = 'reducePca';     % polyZono
    params.R0 = polyZonotope(params.R0);            % convert to polyZono
    params.R = params.R0;                           % for start set
    options.redFactor = 0.0001;                     % zeta_Z (pZ order)
    % zeta_h (timeStep) ... depends on order (chosen in code) and decrFactor
    options.zetaphi = [0.80; 0.75; 0.63];
    options.tensorOrder = 3;                        % fixed
end
options.R.error = zeros(sys.dim,1);                 % for consistency

% options to speed up tensor computation
options.thirdOrderTensorempty = false;
options.isHessianConst = false;
options.hessianCheck = false;

end

function [params,options] = set_nonlinDASys_adaptive(sys,params,options)

options.decrFactor = 0.90;              % adaptation of delta t
options.zetaTlin = 0.0005;              % zeta_T,lin (taylorTerms)
options.zetaTabs = 0.005;               % zeta_T,abs (taylorTerms)
params.R0 = zonotope(params.R0);    % convert to zono
params.R = params.R0;               % for initial set
options.redFactor = 0.0005;         % zeta_Z (zonotope order)
options.zetaK = 0.90;               % zeta_K (tensorOrder)
options.zetaphi = 0.85;             % zeta_Delta (timeStep)
options.R.error = zeros(sys.dim,1);                 % for consistency

% options to speed up tensor computation
options.thirdOrderTensorempty = false;
options.isHessianConst = false;
options.hessianCheck = false;

end

function [params,options] = set_nonlinearSysDT_adaptive(sys,params,options)

if contains(options.alg,'lin')
    options.redFactor = 0.0005; % zeta_Z
    options.tensorOrder = 2;    % just starting value
    options.zetaK = 0.8;        % zeta_K 
elseif contains(options.alg,'poly')
    % options for poly (fixed in adaptive)
    options.polyZono.volApproxMethod = 'interval';           % polyZono
    options.polyZono.maxDepGenOrder = 20;                    % polyZono (unused...)
    options.polyZono.maxPolyZonoRatio = 0.01;                % polyZono
    options.polyZono.restructureTechnique = 'reduceGirard';   % polyZono
    params.R0 = polyZonotope(params.R0);            % convert to polyZono
    params.R = params.R0;                           % for start set
    options.redFactor = 0.0001;                     % zeta_Z (pZ order)
    % zeta_h (timeStep) ... depends on order (chosen in code) and decrFactor
    options.tensorOrder = 3;                        % fixed
end

end

function [params,options] = set_unsatIntervals(sys,params,options)
% returns time intervals, where unsafe sets / safe sets need to be checked
% for computational efficiency, these are stored as matrices and not as
% object of the interval-class

options.savedata.unsafeSet_unsat = cell(length(params.unsafeSet),1);
for i=1:length(params.unsafeSet)
    if ~isempty(params.unsafeSet{i}.time)
        options.savedata.unsafeSet_unsat{i} = ...
            [infimum(params.unsafeSet{i}.time), ...
            supremum(params.unsafeSet{i}.time)];
    else
        % default: unsatisfied over whole time horizon
        options.savedata.unsafeSet_unsat{i} = [params.tStart, params.tFinal];
    end
end

options.savedata.safeSet_unsat = cell(length(params.safeSet),1);
for i=1:length(params.safeSet)
    if ~isempty(params.safeSet{i}.time)
        options.savedata.safeSet_unsat{i} = ...
            [infimum(params.safeSet{i}.time), ...
            supremum(params.safeSet{i}.time)];
    else
        % default: unsatisfied over whole time horizon
        options.savedata.safeSet_unsat{i} = [params.tStart, params.tFinal];
    end
end
    
end

function [params,options] = set_safeUnsafeSetFastInner(sys,params,options)
% if any single safe or unsafe set is just a halfspace perpendicular to an
% axis of the output, we can use a simpler method to check for an
% intersection of the inner-approximation with that set where the
% inner-approximation does not actually have to be computed

% loop over unsafe sets
for i=1:length(params.unsafeSet)
    params.unsafeSet{i}.fastInner = false;
    oneHalfspace = size(params.unsafeSet{i}.set.P.A,1) == 1;
    perpHalfspace = nnz(params.unsafeSet{i}.set.P.A) == 1;
    if oneHalfspace && perpHalfspace
        params.unsafeSet{i}.fastInner = true;
    end
end

% loop over safe sets
for i=1:length(params.safeSet)
    params.safeSet{i}.fastInner = false;
    oneHalfspace = size(params.safeSet{i}.set.P.A,1) == 1;
    perpHalfspace = nnz(params.safeSet{i}.set.P.A) == 1;
    if oneHalfspace && perpHalfspace
        params.safeSet{i}.fastInner = true;
    end
end

end


% Auxiliary Functions: hybridAutomaton and parallelHybridAutomaton

function [params,options] = set_HA_timeStep(sys,params,options)
% timeStep becomes timeStepLoc for each location ... only if non-adaptive

if ~strcmp(options.linAlg,'adaptive')
    
    if ~iscell(options.timeStep)
        locations = sys.location;
        numLoc = length(locations);
        options.timeStepLoc = cell(numLoc,1);
        for i=1:numLoc
            options.timeStepLoc{i} = options.timeStep;
        end
    else
        options.timeStepLoc = options.timeStep;
    end
    
    options = rmfield(options,'timeStep');
    
end

end

function [params,options] = set_HA_inputs(sys,params,options)
% internal value options.Uloc stores the input set for each location

locations = sys.location;
numLoc = length(locations);

if ~iscell(params.U)
    isaInt = false;
    if isa(params.U,'interval')
        isaInt = true; UZon = zonotope(params.U);
    end
    % same input set for each location
    Uloc = cell(numLoc,1);
    for i = 1:numLoc
        if isaInt
            Uloc{i} = UZon;
        else
            Uloc{i} = params.U;
        end
    end
    
else
    % copy input set, convert to zonotope if necessary
    Uloc = cell(numLoc,1);
    for i = 1:numLoc
        if isa(params.U{i},'interval')
            Uloc{i} = zonotope(params.U{i});
        else
            Uloc{i} = params.U{i};
        end
    end
end

params = rmfield(params,'U');
params.Uloc = Uloc;

end

function [params,options] = set_HA_sim_inputs(sys,params,options)
% internal value options.uLoc which the input set for each location

locations = sys.location;
numLoc = length(locations);

if ~iscell(params.u) 
    % same input set for each location    
    uloc = cell(numLoc,1);
    for i = 1:numLoc
        uloc{i} = params.u;
    end
else
    % copy input set for each location
    uloc = params.u;
end

params = rmfield(params,'u');
params.uLoc = uloc;

end

function [params,options] = set_pHA_inputs(sys,params,options)
% internal value options.Uloc stores the input set for each location

numComp = length(sys.components);

if ~iscell(params.U)
    isaInt = false;
    if isa(params.U,'interval')
        isaInt = true; UZon = zonotope(params.U);
    end
    % same input set for each location
    uloc = cell(numComp,1);
    for i = 1:numComp
        numLoc = length(sys.components(i).location);
        uloc{i} = cell(numLoc,1);    
        for j = 1:numLoc
            if isaInt
                uloc{i}{j} = UZon;
            else
                uloc{i}{j} = params.U;
            end
        end
    end
else
    % convert to zonotope 
    for i = 1:numComp
        numLoc = length(sys.components(i).location);
        for j = 1:numLoc
            if isa(params.U{i}{j},'interval')
                uloc{i}{j} = zonotope(params.U{i}{j});
            else
                uloc{i}{j} = params.U{i}{j};
            end
        end
    end    
end

params.Uloc = uloc;

params = rmfield(params,'U');

end

function [params,options] = set_pHA_sim_inputs(sys,params,options)
% internal value options.uloc stores the input set for each location

numComps = length(sys.components);

if ~iscell(params.u) 
    % same input set for each location
    uloc = cell(numComps,1);
    numLoc = length(sys.components(1).location);
    uloc{1} = cell(numLoc,1);
    for i = 1:numLoc
        uloc{1}{i} = params.u;
    end   
    params.uLoc = uloc;
else
    params.uLoc = params.u;
end

params = rmfield(params,'u');

end



%------------- END OF CODE --------------
