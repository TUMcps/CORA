function simRes = simulateRandom(HA,params,varargin)
% simulateRandom - simulates a hybrid automata for random initial points
%                  and random inputs
%
% Syntax:  
%    res = simulateRandom(HA,params,options)
%
% Inputs:
%    HA - hybridAutomaton object
%    params - system parameters
%    options - settings for random simulation
%       .points - nr of simulation runs
%       .fracVert - fraction of initial states starting from vertices
%       .fracInpVert - fraction of input values taken from the 
%                       vertices of the input set
%       .nrConstInp - number of piecewise constant inputs
%
% Outputs:
%    simRes - simResult object storing simulation results
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Niklas Kochdumper
% Written:       03-December-2019 
% Last update:   08-May-2020 (MW, update interface)
% Last revision: ---

%------------- BEGIN CODE --------------

% input argument validation
options = struct();
if nargin == 3 && isstruct(varargin{1})
    options = varargin{1};
end

% options preprocessing
options = validateOptions(HA,mfilename,params,options);

% initialize random inputs
u = cell(size(options.Uloc));

for j = 1:length(u)
    u{j} = randPoint(options.Uloc{j});
end

% determine random points inside the initial set
nrEx = ceil(options.points*options.fracVert);
nrNor = options.points - nrEx;
points = [];
if nrEx > 0
    points = [points, randPoint(options.R0,nrEx,'extreme')]; 
end
if nrNor > 0
    points = [points, randPoint(options.R0,nrNor,'standard')];
end

time = linspace(options.tStart,options.tFinal,options.nrConstInp);
startLoc = options.startLoc;

% initialization
loc = cell(options.points,1);
x = cell(options.points,1);
t = cell(options.points,1);

% simulate the hybrid automaton
for i = 1:options.points

    counter = 1;
    locCur = startLoc;
    
    % loop over all input changes
    for g = 1:length(time)-1
    
        % compute random input
        optsSim.u = u;
        
        if counter < options.nrConstInp * options.fracInpVert 
            optsSim.u{locCur} = randPoint(options.Uloc{locCur},1,'extreme');
        else
            optsSim.u{locCur} = randPoint(options.Uloc{locCur});
        end              
        
        % simulate hybrid automaton
        optsSim.x0 = points(:,i);
        if time(g) ~= 0
            optsSim.tStart = time(g);
        end
        optsSim.tFinal = time(g+1);
        optsSim.startLoc = locCur;
        optsSim.finalLoc = options.finalLoc;
        
        [tTemp,xTemp,locTemp] = simulate(HA,optsSim);
        
        % store results
        [t,x,loc] = aux_concatenateSims(t,x,loc,tTemp,xTemp,locTemp,i);
        % original version:
        % t = [t; tTemp]; x = [x; xTemp]; loc = [loc; locTemp];
        
        % update location and initial point
        points(:,i) = xTemp{end}(end,:)';
        locCur = locTemp{end};
        optsSim = [];
        counter = counter + 1;
    end

end

% for identity resets, we don't need NaN to indicate jumps
[t,x,loc] = aux_identityJumps(t,x,loc);
    
% create object storing the simulation results
simRes = simResult(x,t,loc);

end

% Auxiliary functions -----------------------------------------------------

function [t,x,loc] = aux_concatenateSims(t,x,loc,tTemp,xTemp,locTemp,run)

% state dimension
n = size(xTemp{1},2);

% variables for current part of the trajectory
t_add = []; x_add = []; loc_add = [];

% loop over all parts of the added trajectory to current run
for i=1:length(tTemp)
    % number of steps
    steps = length(tTemp{i});
    % add NaN to indicate potential jump in the trajectory
    t_add = [tTemp{i}; NaN];
    x_add = [xTemp{i}; NaN(1,n)];
    % repeat location to match length of time vecttor
    loc_add = [repmat(locTemp{i},steps,1); NaN];
end

% append to current run
t{run} = [t{run};t_add];
x{run} = [x{run};x_add];
loc{run} = [loc{run};loc_add];

end

function [t,x,loc] = aux_identityJumps(t,x,loc)
% NaN indicates a jump from one location to another, if the reset
% function is the identity matrix, then we can remove the NaN

% loop over all runs
for i=1:length(t)
    % index of kept entries (always remove last entry = NaN)
    idxKeep = [true(length(t{i})-1,1);false];
    
    % loop over all NaNs (except last one, already done above)
    allNaN = find(isnan(t{i}));
    for j=1:length(allNaN)-1
        % explicitly deine meaning of indices
        idxBeforeJump = allNaN(j)-1;
        idxAfterJump = allNaN(j)+1;

        % compare state vector before and after reset
        if all(withinTol(x{i}(idxBeforeJump,:),x{i}(idxAfterJump,:),1e-14))
            % same state vector before and after reset -> remove
            idxKeep(idxBeforeJump:idxBeforeJump+1) = false;
        end
    end

    % remove rows from x, t, loc
    t{i} = t{i}(idxKeep);
    x{i} = x{i}(idxKeep,:);
    loc{i} = loc{i}(idxKeep);

end

end

%------------- END OF CODE --------------