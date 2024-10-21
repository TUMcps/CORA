function R = reachInnerMinkdiff(nlnsys,params,options)
% reachInnerMinkdiff - computes an inner-approximation of the reachable set
%    with the algorithm in [1]
%
% Syntax:
%    R = reachInnerMinkdiff(nlnsys,params,options)
%
% Inputs:
%    nlnsys - nonlinearSys object
%    params - parameters defining the reachability problem
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - object of class reachSet storing the inner-approximation of the 
%        reachable set
%
% References:
%    [1] M. Wetzlinger, A. Kulmburg, and M. Althoff. "Inner approximations
%        of reachable sets for nonlinear systems using the Minkowski
%        difference". IEEE Control Systems Letters, 2024.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSys/reachInner

% Authors:       Mark Wetzlinger
% Written:       17-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% options preprocessing
[params,options] = validateOptions(nlnsys,params,options);

% time period and number of steps
tVec = params.tStart:options.timeStep:params.tFinal;
steps = length(tVec) - 1;

% init struct for reachable sets
timePoint.set = cell(steps+1,1);
timePoint.set{1} = params.R0;
timePoint.time = num2cell(tVec');
timeInt.set = cell(steps,1);
timeInt.time = cellfun(@(x) interval(x(1),x(2)), ...
    num2cell([tVec(1:end-1)',tVec(2:end)'],2),'UniformOutput',false);

% compute derivatives
derivatives(nlnsys,options);

% init step
options.empty = false;
Rtp.set = params.R0;
Rtp.error = zeros(nlnsys.nrOfStates,1);
params.uTrans = params.u;

% check if initial set is a parallelotope (allows for a faster computation)
options.R0ispara = representsa_(params.R0,'interval',1e-10);


% loop over all reachability steps
for k=1:steps
    options.k = k;

    % increment time and log information
    options.t = tVec(k);
    verboseLog(options.verbose,k,options.t,params.tStart,params.tFinal);

    % compute next reachable set
    [Rtp,Rti,options] = linReachInner(nlnsys,Rtp,params,options);
    % exit if inner approximation has become empty
    if options.empty
        break
    end
    
    % compute output set and save in cell structure
    timeInt.set{k} = outputSet(nlnsys,Rti,params,options);
    timePoint.set{k+1} = outputSet(nlnsys,Rtp.set,params,options);
    
end

% construct reachset object
R = reachSet.initReachSet(timePoint,timeInt);

% ------------------------------ END OF CODE ------------------------------
