clc; clear; close all;

%% System dynamics --------------------------------------------------------

dt = 0.1;
A = [1 dt 0 0; 0 1 0 0; 0 0 1 dt; 0 0 0 1];
B = [0.5*dt^2 0; dt 0; 0 0.5*dt^2; 0 dt];
% instantiate discrete-time system
sys = linearSysDT('car',A,B,dt);


%% Model parameters -------------------------------------------------------

params.tStart = 0;
params.tFinal = 1;
params.U = zonotope(interval([-6;-0.2],[6;0.2]));

% assignment to params.R0 later
R0 = zonotope([[10;5;0;0],0.1*diag(ones(sys.dim,1))]);


%% Settings for reachability and simulation

options.zonotopeOrder = 40;
simOpt.points = 20;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 10;

%% Reachability analysis and simulation -----------------------------------

shift = {[0;0;0;0],[25;0;0;0],[50;0;0;0]};
nrCases = length(shift);
R = cell(nrCases,1);
simRes = cell(nrCases,1);

for i=1:nrCases
    % shift initial set
    params.R0 = R0 + shift{i};
    % reachability analysis
    R{i} = reach(sys,params,options);
    % simulation
    simRes{i} = simulateRandom(sys, params, simOpt);
end


%% Visualization ----------------------------------------------------------
% 1 --> s_zeta
% 2 --> v_zeta
% 3 --> s_eta
% 4 --> v_eta

% dims = {[1 2],[3 4]};
dims = {[1 3]};

color_blue = [0, 92, 171] ./ 255;
color_red = [227, 27, 35] ./ 255;

for k = 1:length(dims)

    figure; hold on;
    projDims = dims{k};

    % loop over different initial sets

    for i=1:nrCases
        % plot reachable set
        plot(R{i},projDims,'FaceColor',color_blue,'EdgeColor',color_blue);
    
        % plot simulation
        plot(simRes{i},projDims,'y','LineWidth',2);
    
        % plot initial set
        plot(R0 + shift{i},projDims,'Color',color_red,'LineWidth',2);
    end

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
    xlim([0 80]);
    ylim([-0.4 0.4]);

    grid on;
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperType', 'A4');
    set(gcf, 'PaperPositionMode','manual');
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    orient landscape;
end