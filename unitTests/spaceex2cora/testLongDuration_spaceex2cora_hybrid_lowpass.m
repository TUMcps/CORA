function res = testLongDuration_spaceex2cora_hybrid_lowpass()
% testLongDuration_spaceex2cora_hybrid_lowpass - example for hybrid dynamics;
%    in addition to the converted spaceex model to a flat HA,
%    the same model, but converted to a parallel HA,
%    as well as the original lowpass_parallel.m-file
%    are simulated and subsequently compared to one another
%    the simulation results have to be within a given numerical precision
%
% Syntax:  
%    testLongDuration_spaceex2cora_hybrid_lowpass
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
% 
% Author:       Mark Wetzlinger
% Written:      7-December-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

res = false;

%% automaton #1: original file

PHA = lowpassFilter();

% parameter
params.tFinal = 0.4;
params.x0 = [0;0;0;0]; 
params.startLoc = [1;3]; 

% simulation
[~,simRes{1}] = simulate(PHA,params);


%% automaton #2: converted parallel HA

spaceex2cora('lowpass_parallel.xml',1,[],'lowpass_parallel');
HA_SX = lowpass_parallel();

% simulation 
[~,simRes{2}] = simulate(HA_SX, params);


%% automaton #3: converted flat HA

spaceex2cora('lowpass_parallel.xml',0,[],'lowpass_flat');
HA_SXflat = lowpass_flat();

% simulation
params.startLoc = 3;
[~,simRes{3}] = simulate(HA_SXflat, params);


%% compare simulation results
tol = 1e-6;

for i = 1:length(simRes{1})
   if ~all(size(simRes{1}{i}) == size(simRes{2}{i})) || ...
      max(max(abs(simRes{1}{i} - simRes{2}{i}))) > tol 
        error(['Results for original and converted parallel hybrid ', ...
               'automaton are different!']);
   end
   if ~all(size(simRes{1}{i}) == size(simRes{3}{i})) || ... 
      max(max(abs(simRes{1}{i} - simRes{3}{i}))) > tol
        error(['Results for original and converted flat hybrid ', ...
               'automaton are different!']);
   end
end

res = true;

end

%------------- END OF CODE --------------
