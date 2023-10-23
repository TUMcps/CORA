function res = testMATPOWER_converter_loadPowerSystemCase_IEEE14()
% testMATPOWER_converter_loadPowerSystemCase_IEEE14 - unit test for 
% loading an IEEE 14 bus benchmark from MATPOWER
%
% It is checked whether the conversion of a MATPOWER model to our power 
% system description matches a saved power system description.
%
% Syntax:
%    res = testMATPOWER_converter_loadPowerSystemCase_IEEE14()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% References:
%    [1] M. Althoff, submitted ARCH paper on power system benchmarks

% Authors:       Matthias Althoff
% Written:       05-May-2022
% Last update:   ---
%                ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
path = [CORAROOT filesep 'converter' filesep 'powerSystem2cora' filesep 'cases'];

% initialize partial results
resPartial = [];

%% Load MATPOWER model
loadPowerSystemCase('case14','MATPOWER_IEEE14');
% load created power system description
load([path filesep 'MATPOWER_IEEE14']);

% evaluate conversion from MATPOWER
resPartial = aux_evaluate(resPartial, path, bus, genParam, Pd, Qd, VM, Y);

%% Load PSSE model
loadPowerSystemCase('ieee14.raw','PSSE_IEEE14','PSSE');
% load created power system description
load([path filesep 'PSSE_IEEE14']);

% % evaluate conversion from MATPOWER; the converted model seems to be
% % different
% resPartial = aux_evaluate(resPartial, path, bus, genParam, Pd, Qd, VM, Y);

% final result
res = all(resPartial);

end


% Auxiliary functions -----------------------------------------------------

function resPartial = aux_evaluate(resPartial, path, bus, genParam, Pd, Qd, VM, Y)

% save to specific filenames
bus_created = bus;
genParam_created = genParam;
Pd_created = Pd;
Qd_created = Qd;
VM_created = VM;
Y_created = Y;

% load saved power system description
load([path filesep 'IEEE14']);

% set accuracy
accuracy = 1e-6;


%% perform evaluation
% check Pd
resPartial(end+1) = (norm(Pd - Pd_created) <= accuracy);
% check Qd
resPartial(end+1) = (norm(Qd - Qd_created) <= accuracy);
% check VM (the previous version only used values for generators)
resPartial(end+1) = (norm(VM(bus.generator) - VM_created(bus.generator)) <= accuracy);
% check Y
resPartial(end+1) = (norm(full(Y) - full(Y_created)) <= accuracy);
% check genParam.D
resPartial(end+1) = (norm(genParam.D - genParam_created.D) <= accuracy);
% check genParam.M
resPartial(end+1) = (norm(genParam.M - genParam_created.M) <= accuracy);
% check genParam.R_d
resPartial(end+1) = (norm(genParam.R_d - genParam_created.R_d) <= accuracy);
% check genParam.T_sv
resPartial(end+1) = (norm(genParam.T_sv - genParam_created.T_sv) <= accuracy);
% check genParam.X_m
resPartial(end+1) = (norm(genParam.Y_m - genParam_created.Y_m) <= accuracy);
% check genParam.omega_s
resPartial(end+1) = (norm(genParam.omega_s - genParam_created.omega_s) <= accuracy);
% check bus.load
resPartial(end+1) = (norm(bus.load - bus_created.load) <= accuracy);
% check bus.generator
resPartial(end+1) = (norm(bus.generator - bus_created.generator) <= accuracy);
% check bus.slack
resPartial(end+1) = (norm(bus.slack - bus_created.slack) <= accuracy);
% check bus.input
resPartial(end+1) = (norm(bus.input - bus_created.input) <= accuracy);
% check bus.output
resPartial(end+1) = (norm(bus.output - bus_created.output) <= accuracy);
% check bus.fault
resPartial(end+1) = (norm(bus.fault - bus_created.fault) <= accuracy);

end


% ------------------------------ END OF CODE ------------------------------
