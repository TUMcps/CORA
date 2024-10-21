function [R, res, Verror] = reach(nlnsysDT,params,options,varargin)
% reach - computes the reachable sets of the discrete time system
%
% Syntax:
%    [R,res,Verror] = reach(nlnsysDT,params,options)
%    [R,res,Verror] = reach(nlnsysDT,params,options,spec)
%
% Inputs:
%    nlnsysDT - nonlinearSysDT object
%    params - parameter defining the reachability problem
%    options - options for the computation of the reachable set
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res - true/false whether specifications are satisfied
%    Verror - linearization error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       21-August-2012
% Last update:   29-January-2018
%                19-November-2022 (MW, integrate output equation)
%                10-May-2023 (LL, integrate uTrans in U)
%                06-November-2023 (LL, add Verror as output)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
spec = setDefaultValues({[]},varargin);

% options preprocessing
[params,options] = validateOptions(nlnsysDT,params,options);

% compute symbolic derivatives
derivatives(nlnsysDT,options);

% initialize cell array that stores the reachable sets
t = params.tStart:nlnsysDT.dt:params.tFinal;
steps = length(t)-1;
timePoint.set = cell(steps+1,1);
Verror.L_y = cell(steps+1,1);
Verror.L_x = cell(steps+1,1);

% add constant input
if isfield(params,'uTrans')
    params.U = params.U + params.uTrans;
    params.uTrans = 0;
end
U0 = params.U;

% add input for time 1
if isfield(params,'uTransVec')
    params.U = U0 + params.uTransVec(:,1);
end  

% compute output for time 1
[timePoint.set{1}, Verror.L_y{1}] = outputSet(nlnsysDT,params.R0,params,options);
Rnext = params.R0;
Verror.L_x{1} = zonotope(zeros(size(Rnext.c)));

% loop over all reachablity steps
for i = 1:steps

    options.i = i;
    % compute next reachable set
    [Rnext,options,Verror.L_x{i+1}] = linReach(nlnsysDT,Rnext,params,options);

    % add input for time i if a trajectory should be tracked
    if isfield(params,'uTransVec')
        params.U = U0 + params.uTransVec(:,i+1);
    end  
    % compute output set
    [timePoint.set{i+1}, Verror.L_y{i+1}] = outputSet(nlnsysDT,Rnext,params,options);

    % log information
    verboseLog(options.verbose,i,t(i),params.tStart,params.tFinal);
    
    % check specification
    if ~isempty(spec)
       if ~check(spec,timePoint.set{i+1},interval(t(i+1)))
           timePoint.set = timePoint(1:i+1);
           timePoint.time = num2cell(t(1:i+1)');
           R = reachSet(timePoint);
           res = false;
           return;
       end
    end
end

% create reachable set object
timePoint.time = num2cell(t(1:end)');
R = reachSet(timePoint);

% log information
verboseLog(options.verbose,length(t),t(end),params.tStart,params.tFinal);

% ------------------------------ END OF CODE ------------------------------
