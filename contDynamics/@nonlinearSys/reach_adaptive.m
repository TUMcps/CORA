function [timeInt,timePoint,res,tVec,options] = reach_adaptive(obj,params,options)
% reach_adaptive - computes the reachable continuous set
%
% Syntax:
%    [timeInt,timePoint,res,tVec,options] = reach_adaptive(obj,params,options)
%
% Inputs:
%    obj - continuous system object
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    timeInt - cell-array of time-interval solutions
%    timePoint - cell-array of time-point solutions
%    res - satisfaction / violation of specifications
%    tVec - vector of time steps
%    options - options for the computation of reachable sets (param tracking)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       27-May-2020
% Last update:   25-February-2021 (merge to master)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize cell-arrays that store the reachable set
timeInt.set = {};
timeInt.time = {};
timePoint.set = {options.R0};
timePoint.time = {options.tStart};
res = 0;
tVec = 0;

% remove 'adaptive' from alg (just for tensor computation)
if contains(options.alg,'lin')
    options.alg = 'lin';
elseif contains(options.alg,'poly')
    options.alg = 'poly';
end

% iteration counter and time for main loop
options.i = 1;
options.t = options.tStart;
% set linearization errors (separate values for Delta and optimal Delta t)
options.error_adm_horizon = zeros(obj.dim,1);
options.error_adm_Deltatopt = zeros(obj.dim,1);
% init abortion flag
abortAnalysis = false;

% just for safety
if ~isfield(options,'linReach_new')
    options.linReach_new = true;
end

% timePoint.set{1,1} = options.R0;
% timePoint.time{1,1} = options.t;

% MAIN LOOP
while options.tFinal - options.t > 1e-12 && ~abortAnalysis
    
    % log information
    verboseLog(options.i,options.t,options);

    % reduction of R via restructuring (only poly)
    if isa(options.R,'polyZonotope')
        ratio = approxVolumeRatio(options.R,options.polyZono.volApproxMethod);
        if ratio > options.polyZono.maxPolyZonoRatio
            options.R = restructure(options.R,...
                options.polyZono.restructureTechnique,options.polyZono.maxDepGenOrder);
        end
    end
    
    % propagation of reachable set
    if options.linReach_new
        [Rnext.ti,Rnext.tp,options] = linReach_adaptive(obj,options,options.R);
    else
        [Rnext.ti,Rnext.tp,options] = linReach_adaptive_old(obj,options,options.R);
    end
    
    % reduction for next step
    Rnext.ti = reduce(Rnext.ti,'adaptive',options.redFactor*5); % not reused
    Rnext.tp = reduce(Rnext.tp,'adaptive',options.redFactor);
    % try: additional reduction for poly
    if isa(Rnext.tp,'polyZonotope')
        Rnext.tp = aux_reduceOnlyDep(Rnext.tp,options.polyZono.maxDepGenOrder);
    end
    % track zonotope orders
    if isa(Rnext.tp,'polyZonotope')
        options.zonordersRtp(options.i,1) = size(Rnext.tp.G,2) / obj.dim;
        options.zonordersRtp(options.i,2) = size(Rnext.tp.GI,2) / obj.dim;
    else
        options.zonordersRtp(options.i,1) = size(generators(Rnext.tp),2) / obj.dim;
    end
    % track abstraction order
    options.kappa(options.i,1) = options.tensorOrder;
    

    % save to output variables
    tVec(options.i,1) = options.timeStep;
    % save reachable set in cell structure
    timeInt.set{options.i,1} = Rnext.ti; 
    timeInt.time{options.i,1} = interval(options.t,options.t+tVec(options.i));
    timePoint.set{options.i+1,1} = Rnext.tp;
    timePoint.time{options.i+1,1} = options.t+tVec(options.i);
    
    % increment time
    options.t = options.t + options.timeStep;
    
    % update iteration counter
    options.i = options.i + 1;
    
    % start set for next step (since always initReach called)
    options.R = Rnext.tp;
    
    % check for timeStep -> 0
    abortAnalysis = aux_checkForAbortion(tVec,options.t,options.tFinal);
end

% log information
verboseLog(options.i,options.t,options);

end


% Auxiliary functions -----------------------------------------------------

function abortAnalysis = aux_checkForAbortion(tVec,currt,tFinal)
% check last N steps of time step vector: if those time steps are too small,
% we expect this to continue so that the analysis will not terminate

% init flag
abortAnalysis = false;

% remaining time
remTime = tFinal - currt;
% number of previous steps considered for criterion
N = 10;
% total steps until now
k = length(tVec);

% criterion: if sum of last N steps is smaller than a certain fraction of
%            the remaining time, abort the analysis
lastNsteps = sum(tVec(end-min(N,k)+1:end));
if remTime / lastNsteps > 1e9
    abortAnalysis = true;
    warning(sprintf(['The analysis is aborted because the time step size converges to 0.\n'...
        '         The reachable sets until t = ' num2str(currt) ' are returned.']));
end

end

function Rnew = aux_reduceOnlyDep(R,order)

[n, Gsize] = size(R.G);
if Gsize / n < order + 1
    Rnew = R;
    return;
end

% order generators by length
h = vecnorm(R.G,2);

% determine the smallest generators (= generators that are removed)
[~,ind] = sort(h,'descend');
ind = ind(order*n+1:end);

% construct a zonotope from the generators that are removed
Gred = R.G(:,ind);
Ered = R.E(:,ind);
pZred = polyZonotope(zeros(n,1),Gred,[],Ered);

% zonotope over-approximation
zono = zonotope(pZred);
zono = zonotope(zono.c,diag(sum(abs(zono.G),2)));

% remove the generators that got reduced from the generator matrices
Grem = R.G;
Grem(:,ind) = [];
Erem = R.E;
Erem(:,ind) = [];

% add the reduced generators as new independent generators 
newc = R.c + zono.c;
GInew = [R.GI, zono.G];

% instantiate new R
Rnew = polyZonotope(newc,Grem,GInew,Erem);

end

% ------------------------------ END OF CODE ------------------------------
