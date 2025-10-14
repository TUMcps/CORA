function [res, R, traj] = verify(obj, spec, params, options, varargin)
% verify - tries to verify the specification with the given params/options
%    1. simulation random rusn
%    2. check violations in simulations
%       2.1 if violated: set initial set to starting point
%       2.2 no violation: keep initial set
%    3. compute reachable set
%
% Syntax:
%    [res, R, traj] = verify(obj, spec, params, options)
%    res = verify(obj, spec, params, options)
%
% Inputs:
%    obj - neurNetContrSys object
%    spec - object of class specification
%    params - parameter defining the reachability problem
%    options - options for the computation of reachable sets and network evaluation
%    verbose - verbose output
%
% Outputs:
%    res - 'VIOLATED', 'VERIFIED', 'UNKNOWN'
%    R - object of class reachSet storing the computed reachable set
%    traj - trajectory object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neurNetContrSys

% Authors:       Tobias Ladner
% Written:       22-March-2023
% Last update:   20-July-2023 (TL, check early stop)
%                02-May-2024 (TL, bug fix violating runs)
%                23-April-2025 (TL, added recursive splitting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(4,6);
[verbose,splitR0] = setDefaultValues({false,0}, varargin);
res = [];

inputArgsCheck({ ...
    {obj, 'att', 'neurNetContrSys'}; ...
    {spec, 'att', 'specification'}; ...
    {params, 'att', 'struct'}; ...
    {options, 'att', 'struct'}; ...
    {verbose, 'att', 'logical', 'scalar'}; ...
    {splitR0, 'att', 'numeric', 'scalar'}; ...
})

% Simulation --------------------------------------------------------------

timerVal = tic;
traj = simulateRandom(obj, params);
tSim = toc(timerVal);
if verbose
    disp(['Time to compute random simulations: ', num2str(tSim)]);
end

% Check Violation ---------------------------------------------------------

timerVal = tic;
[isNotVio, ~, indObj] = check(spec, traj);
isVio = ~isNotVio; % spec is safeSet
tVio = toc(timerVal);
if verbose
    disp(['Time to check violation in simulations: ', num2str(tVio)]);
end

if isVio
    % only continue with violating run
    traj = trajectory([], ...
        traj(indObj{1}).x(:,:,indObj{2}), [], ...
        traj(indObj{1}).t);
    params.R0 = polyZonotope(traj.x(:,1,1));
    
end

% Reachability Analysis ---------------------------------------------------

isRchecked = false;
timerVal = tic;
if isVio
    % do not include spec while verifying violating run
    R = reach(obj, params, options);
else
    % verify entire initial set; stop early in case of violation
    R = reach(obj, params, options, spec);

    % check if initial set should be split
    isVeri = (R(end).timePoint.time{end} == params.tFinal && check(spec, R));
    if ~isVeri && splitR0 > 0
        % reachability analysis has not finished
        % split recursively
        if verbose
            disp('Unable to verify it directly. Trying recursive splitting..')
        end
        [res, R] = aux_VerifyRecursively(obj,spec,params,options,verbose,splitR0);
        isVeri = strcmp(res,'VERIFIED');
    end

    % specs already checked
    isRchecked = true;
end

% stop time
tComp = toc(timerVal);
if verbose
    disp(['Time to compute reachable set: ', num2str(tComp)]);
end

% Verification ------------------------------------------------------------

if ~isRchecked
    timerVal = tic;
    isVeri = check(spec, R);
    tVeri = toc(timerVal);
    if verbose
        disp(['Time to check verification: ', num2str(tVeri)]);
    end
else
    if verbose
        disp('Reachable set already verified during computation.')
    end
    tVeri = 0;
end

% Finish ------------------------------------------------------------------

if isempty(res)
    if isVio && ~isVeri
        res = 'FALSIFIED';
    elseif isVeri
        res = 'VERIFIED';
    else
        res = 'UNKNOWN';
    end
end

if verbose
    tTotal = tSim+tVio+tComp+tVeri;
    disp(['Total Time: ', num2str(tTotal)]);
end

end


% Auxiliary functions -----------------------------------------------------

function [res, R] = aux_VerifyRecursively(obj,spec,params,options,verbose,splitR0)
% verifies the system recursively through input set splitting

% init
table = CORAtable("single",{'Time [s]','Iteration','#Processed','#Queued','#Verified'},{'.1f','i','i','i','i'});
table.printHeader();
listR0 = {params.R0};
numVerified = 0;
numProcessed = 0;
% use (initial) initial set as first set to properly plot it later
timePoint.set = {params.R0,params.R0};
timePoint.time = {0,0};
timeInterval.set = {params.R0};
timeInterval.time = {interval(0)};
R = reachSet(timePoint,timeInterval);

% iterate and split
timerVal = tic;
for i=1:splitR0
    % print content row
    table.printContentRow({toc(timerVal),i,numProcessed,numel(listR0),numVerified});

    % init
    listR0split = cell(1,numel(listR0)*2);
    counter = 0;

    % go through all initial sets
    for j=1:numel(listR0)
        % call reach with current set
        params.R0 = listR0{j};

        % only call reach if i>1 (as already computed)
        if i>1
            R_j = reach(obj, params, options, spec);
        end

        % check specification
        if i>1 && (R_j(end).timePoint.time{end} == params.tFinal) && check(spec, R_j)
            % verified; add to computed reachable set
            numVerified = numVerified + 1;
            R = [R;R_j];
        else
            % split along longest dimension and save for next iteration
            R0 = params.R0;
            [~,n] = max(rad(interval(R0)).*mean(abs(obj.nn.calcSensitivity(center(R0))),1)');
            listR0split(counter+1:counter+2) = split(R0,n);
            counter = counter+2;

            % save reachable set if end of loop is reached
            if i==splitR0
                R = [R;R_j];
            end
        end
    end
    numProcessed = numProcessed + numel(listR0);
    
    % update list
    listR0 = listR0split(1:counter);
    if isempty(listR0)
        break;
    end
end

% finish table
table.printContentRow({toc(timerVal),splitR0+1,numProcessed,numel(listR0),numVerified});
table.printFooter();

% check if all are verified
if isempty(listR0)
    res = 'VERIFIED';
else
    res = 'UNKNOWN';
end

end

% ------------------------------ END OF CODE ------------------------------
