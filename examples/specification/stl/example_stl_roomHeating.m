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

% Warning
disp('Warning: This example typically runs for many hours. Please comment the following lines to run it.');
completed = false;
return

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

tic
R = reach(HA, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Verification ------------------------------------------------------------

T = stl('T',2);

phi{1} = finally(globally(T(1) > 19.8 & T(1) < 20.0 & T(2) > 19.9 ...
    & T(2) < 21.1, interval(0,9)), interval(0,2.2));

res = true;


% 1. Three-valued Signals

tic

for i = 1:length(phi)
    res = res && modelChecking(R,phi{i},'signals');

    if ~res
        throw(CORAerror('CORA:testFailed'));
    end
end

tComp = toc;
disp(['computation time with signals: ',num2str(tComp)]);


% 2. Reachset Temporal Logic
% takes a very long time

tic

for i = 1:length(phi)
    res = res && modelChecking(R,phi{i},'rtl');

    if ~res
        throw(CORAerror('CORA:testFailed'));
    end
end

tComp = toc;
disp(['computation time with rtl: ',num2str(tComp)]);


% 3. Sampled Time
% takes a very long time

tic

for i = 1:length(phi)
    res = res && modelChecking(R,phi{i},'sampledTime');

    if ~res
        throw(CORAerror('CORA:testFailed'));
    end
end

tComp = toc;
disp(['computation time with sampledTime: ',num2str(tComp)]);


end

% ------------------------------ END OF CODE ------------------------------
