function [timeInt,timePoint,res,tVec,options] = reach_adaptive(obj,params,options)
% reach_adaptive - computes the reachable continuous set
%
% Syntax:  
%    [R,res,tVec,options] = reach_adaptive(obj,options)
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

% Author:        Mark Wetzlinger
% Written:       27-May-2020
% Last update:   25-February-2021 (merge to master)
% Last revision: ---

%------------- BEGIN CODE --------------

% initialize cell-arrays that store the reachable set
timeInt.set = {};
timeInt.time = {};
timePoint.set = {};
timePoint.time = {};
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

% MAIN LOOP
while options.tFinal - options.t > 1e-12
    
    % log information
    verboseLog(options.i,options.t,options);


    % reduction of R via restructuring (only poly)
    if isa(options.R,'polyZonotope')
        ratio = approxVolumeRatio(options.R,options.volApproxMethod);
        if ratio > options.maxPolyZonoRatio
            options.R = restructure(options.R,...
                options.restructureTechnique,options.maxDepGenOrder);
        end
    end
    
    [Rnext.ti,Rnext.tp,options] = linReach_adaptive(obj,options,options.R);
    
    % reduction for next step
    Rnext.ti = reduce(Rnext.ti,'adaptive',options.redFactor*5); % not reused
    Rnext.tp = reduce(Rnext.tp,'adaptive',options.redFactor);
    % track zonotope orders
    if isa(Rnext.tp,'polyZonotope')
        options.zonordersRtp(options.i,1) = size(Rnext.tp.G,2) / obj.dim;
        options.zonordersRtp(options.i,2) = size(Rnext.tp.Grest,2) / obj.dim;
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
    timePoint.set{options.i,1} = Rnext.tp;
    timePoint.time{options.i,1} = options.t+tVec(options.i);
    
    % increment time
    options.t = options.t + options.timeStep;
    
    % update iteration counter
    options.i = options.i + 1;
    
    % start set for next step (since always initReach called)
    options.R = Rnext.tp;
end

% log information
verboseLog(options.i,options.t,options);

end

%------------- END OF CODE --------------
