function res = example_stl_roomHeating
% example_stl_roomHeating - example of signal temporal logic
% checking of the room heating model
%
% Syntax:
%    res = example_stl_roomHeating()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean

% Authors:       Benedikt Seidl
% Written:       14-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------

HA = roomHeating();


% Parameter ---------------------------------------------------------------

params.tFinal = 12;
params.startLoc = 1;
params.R0 = zonotope([20.5;20.5],diag([0.1,0.1]));
params.U = zonotope(4,0.01);


% Reachability Settings ---------------------------------------------------

% settings for continuous reachability
options.taylorTerms = 5;
options.zonotopeOrder = 10;
options.timeStep = 0.01;

% settings for hybrid systems
options.enclose = {'box','pca'};
options.guardIntersect = 'zonoGirard';

% Reachability Analysis ---------------------------------------------------

timerVal = tic;
R = reach(HA, params, options);
tComp = toc(timerVal);
disp(['computation time of reachable set: ',num2str(tComp)]);


% Verification ------------------------------------------------------------

T = stl('T',2);

phi{1} = finally(globally(T(1) > 19.8 & T(1) < 20.0 & T(2) > 19.9 ...
    & T(2) < 21.1, interval(0,9)), interval(0,2.2));

res = true;


% 1. Three-valued Signals

timerVal = tic;

for i = 1:length(phi)
    assertLoop(modelChecking(R,phi{i},'signals'),i);
end

tComp = toc(timerVal);
disp(['computation time with signals: ',num2str(tComp)]);

% 2. Incremental with four-valued signals

timerVal = tic;

for i = 1:length(phi)
    assertLoop(modelChecking(R,phi{i},'incremental','propFreq',100,'verbose',true),i);
end

tComp = toc(timerVal);
disp(['computation time with incremental: ',num2str(tComp)]);

% Warning
disp('Warning: The next part of this example typically runs for many hours. Please comment the following lines to run it.');
completed = false;
return

% 3. Reachset Temporal Logic
% takes a very long time

timerVal = tic;

for i = 1:length(phi)
    assertLoop(modelChecking(R,phi{i},'rtl'),i);
end

tComp = toc(timerVal);
disp(['computation time with rtl: ',num2str(tComp)]);


% 4. Sampled Time
% takes a very long time

timerVal = tic;

for i = 1:length(phi)
    assertLoop(modelChecking(R,phi{i},'sampledTime'),i);
end

tComp = toc(timerVal);
disp(['computation time with sampledTime: ',num2str(tComp)]);


end

% ------------------------------ END OF CODE ------------------------------
