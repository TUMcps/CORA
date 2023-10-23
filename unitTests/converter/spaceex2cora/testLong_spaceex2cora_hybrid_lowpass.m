function res = testLong_spaceex2cora_hybrid_lowpass()
% testLong_spaceex2cora_hybrid_lowpass - example for hybrid dynamics;
%    in addition to the converted spaceex model to a flat HA,
%    the same model, but converted to a parallel HA,
%    as well as the original lowpass_parallel.m-file
%    are simulated and subsequently compared to one another
%    the simulation results have to be within a given numerical precision
%
% Syntax:
%    res = testLong_spaceex2cora_hybrid_lowpass
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Mark Wetzlinger
% Written:       07-December-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% automaton #1: original file

pHA = lowpassFilter();

% parameter
params.tFinal = 0.4;
params.x0 = [0;0;0;0]; 
params.startLoc = [1;3]; 

% simulation
[~,simRes{1}] = simulate(pHA,params);


%% automaton #2: converted parallel HA

spaceex2cora('lowpass_parallel.xml',true,[],'lowpass_parallel');
pHA_SX = lowpass_parallel();

% simulation 
[~,simRes{2}] = simulate(pHA_SX, params);


%% automaton #3: converted flat HA

spaceex2cora('lowpass_parallel.xml',false,[],'lowpass_flat');
HA_SX = lowpass_flat();

% simulation
params.startLoc = 3;
[~,simRes{3}] = simulate(HA_SX, params);


%% visualization

% figure;
% projDims = {[1,2],[3,4]};
% 
% for p=1:length(projDims)
%     subplot(1,2,p); hold on; box on;
%     for i=1:length(simRes{1})
%         plot(simRes{1}{i}(:,projDims{p}(1)),simRes{1}{i}(:,projDims{p}(2)),...
%             'Color',colorblind('b'));
%         plot(simRes{2}{i}(:,projDims{p}(1)),simRes{2}{i}(:,projDims{p}(2)),...
%             'Color',colorblind('r'),'LineStyle','--');
%         plot(simRes{3}{i}(:,projDims{p}(1)),simRes{3}{i}(:,projDims{p}(2)),...
%             'Color',colorblind('y'),'LineStyle','-.');
%     end
% end


%% compare simulation results
tol = 1e-6;

for i = 1:length(simRes{1})
    % compare results for original and converted parallel hybrid automaton
    if ~all(size(simRes{1}{i}) == size(simRes{2}{i})) || ...
            ~all(all(withinTol(simRes{1}{i},simRes{2}{i},tol)))
        throw(CORAerror('CORA:testFailed'));
    end
    % compare results for original and converted flat hybrid automaton
    if ~all(size(simRes{1}{i}) == size(simRes{3}{i})) || ... 
            ~all(all(withinTol(simRes{1}{i},simRes{2}{i},tol)))
        throw(CORAerror('CORA:testFailed'));
    end
end

res = true;

end

% ------------------------------ END OF CODE ------------------------------
