function verboseLog(verbose,step,t,tStart,tFinal)
% verboseLog - standardized console output if options.verbose = true,
%    only used in reach-functions
%
% Syntax:
%    verboseLog(step,t,options)
%
% Inputs:
%    verbose - true/false
%    step - current step
%    t - start time of current step
%    tStart - start time of time horizon
%    tFinal - end time of time horizon
%
% Outputs:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       05-March-2021
% Last update:   27-June-2022 (MW, add handling for hybrid systems)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~verbose
    return
end

% check if this is called from hybrid system
st = dbstack("-completenames");
for i=1:length(st)
    if contains(st(i).file,['@location' filesep 'reach'])
        % no log in contDynamics evaluation for hybrid systems
        return;
    end
end

% start message
if abs(t - tStart) < 1e-12
    disp(newline + "Start analysis...");
    disp("- step 1: " + t);
    return;
end

% end message
if abs(t - tFinal) < 1e-12
    disp("...time horizon reached, analysis finished." + newline);
    return;
end

% every ... steps, information is logged
cycle = 10;
% shift by one to obtain clean start times
if mod(step-1,cycle) == 0
    disp("- step " + step + ": " + t);
end

% ------------------------------ END OF CODE ------------------------------
