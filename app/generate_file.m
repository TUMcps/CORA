function [id, err] = generate_file(file, handles, currentOption,func_name)
% the err variable is needed to stop the code from running in case there
% was any error in the assigned values of the parameters
err = 0;
id = fopen(file, 'wt');

% in case there was any problem with opening the file
if id == -1
    return
end

disp('CORA APP:')
disp('- Generating script file..')
fprintf('  %s\n',file)

% defining some parameters that will be used later depending on the system
switch currentOption
    case 'Linear System'
        
        % name of the function
        func         = 'Linear';
        
        % the dimensions and nrOfInputs of the system are needed to make sure 
        % the other parameters are defined correctly 
        dimensions    = handles.dimensions_Linear;
        nrOfInputs    = handles.nrOfInputs_Linear;
        Sys           = 'linSys';
        
        % settings parameters
        tstart        = get(handles.txtTstart_Linear, 'String');
        tfinal        = get(handles.txtTfinal_Linear, 'String');
        time_step     = get(handles.txtTimeStep_Linear, 'String');
        
        % algorithim settings
        contents      = cellstr(get(handles.popLA_Linear, 'String'));
        algorithm     = contents{get(handles.popLA_Linear, 'Value')};
        OptError      = get(handles.txtOptionError, 'String');
        zOrder        = get(handles.txtZOrder_Linear, 'String');
        taylor_terms  = get(handles.txtTaylor_Linear, 'String');
        contents      = cellstr(get(handles.popRT_Linear, 'String'));
        RT            = contents{get(handles.popRT_Linear, 'Value')};
        
        % initial and input sets
        R0_handle     = 'R0';
        U_handle      = 'U';
        uTrans_handle = 'uTrans';
        
        % simulation settings
        Simulation    = get(handles.rbSim_Linear, 'Value');
        x0            = get(handles.txtX0_Sim_Linear, 'String');
        u             = get(handles.txtU_Sim_Linear, 'String');
        
        % random simulation settings
        Random_Simulation         = get(handles.rbRSim_Linear, 'Value');
        simOpt_points_random      = get(handles.txtNoP_RSim_Linear, 'String');
        simOpt_fracVert_random    = get(handles.txtFV_RSim_Linear, 'String');
        simOpt_fracInpVert_random = get(handles.txtFIV_RSim_Linear, 'String');
        simOpt_nrConstInp_random  = get(handles.txtNCI_RSim_Linear, 'String');
        
        % RRT simulation settings
        Simulate_RRT          = get(handles.rbSimRRT_Linear, 'Value');
        simOpt_points_RRT     = get(handles.txtNoP_RSim_Linear, 'String');
        Yes                   = get(handles.rbYesEPS_SimRRT_Linear, 'Value');
        No                    = get(handles.rbNoEPS_SimRRT_Linear, 'Value');
        simOpt_stretchFac_RRT = get(handles.txtSF_SimRRT_Linear, 'String');
        
        % plotting settings
        dims = get(handles.listPlots_Linear, 'String');
        
        % reachable set plot settings
        reach_plot                  = get(handles.checkReachPlot_Linear, 'Value');
        color_contents_reach        = cellstr(get(handles.popColorReach_Linear, 'String'));
        color_reach                 = color_contents_reach{get(handles.popColorReach_Linear, 'Value')};
        edge_color_contents_reach   = cellstr(get(handles.popEdgeColorReach_Linear, 'String'));
        edge_color_reach            = edge_color_contents_reach{get(handles.popEdgeColorReach_Linear, 'Value')};
        
        % initial set plot settings
        initial_plot                = get(handles.checkInitialPlot_Linear, 'Value');
        color_contents_initial      = cellstr(get(handles.popColorInitial_Linear, 'String'));
        color_initial               = color_contents_initial{get(handles.popColorInitial_Linear, 'Value')};
        edge_color_contents_initial = cellstr(get(handles.popEdgeColorInitial_Linear, 'String'));
        edge_color_initial          = edge_color_contents_initial{get(handles.popEdgeColorInitial_Linear, 'Value')};

        % simulation plot settings
        simulation_plot             = get(handles.checkSimulationPlot_Linear, 'Value')  && strcmp(get(handles.checkSimulationPlot_Linear, 'Enable'),'on');
        color_contents_simulation   = cellstr(get(handles.popColorSimulation_Linear, 'String'));
        color_simulation            = color_contents_simulation{get(handles.popColorSimulation_Linear, 'Value')};
        line_contents_simulation    = cellstr(get(handles.popLineSimulation_Linear, 'String'));
        line_simulation             = line_contents_simulation{get(handles.popLineSimulation_Linear, 'Value')};   
        
        % auxiliary function
        SpaceEx = get(handles.rbLoadSModel_Linear, 'Value');
        
    case 'Nonlinear System'
        
        % name of the function
        func         = 'Nonlinear';

        % settings parameters
        tstart        = get(handles.txtTstart_Nonlinear, 'String');
        tfinal        = get(handles.txtTfinal_Nonlinear, 'String');
        time_step     = get(handles.txtTimeStep_Nonlinear, 'String');
        
        % initial and input sets 
        R0_handle     = 'R0_Nonlinear';
        U_handle      = 'U_Nonlinear';
        uTrans_handle = 'uTrans_Nonlinear';
        Sys           = 'nonlinSys';
        
        % algorithim settings
        zOrder        = get(handles.txtZOrder_Nonlinear, 'String');
        taylor_terms  = get(handles.txtTaylor_Nonlinear, 'String');
        contents      = cellstr(get(handles.popRT_Nonlinear, 'String'));
        RT            = contents{get(handles.popRT_Nonlinear, 'Value')};
        
        % conservative linearization
        tensor_order_lin        = get(handles.txtConLinTO_Nonlinear, 'String');
        intermediate_order_lin  = get(handles.txtIO_Nonlinear_lin, 'String');
        error_order_lin         = get(handles.txtEO_Nonlinear_lin, 'String');
        
        % conservative polynomilization
        tensor_order_poly       = get(handles.txtConPolTO_Nonlinear, 'String');
        intermediate_order_poly = get(handles.txtIO_Nonlinear, 'String');
        error_order_poly        = get(handles.txtEO_Nonlinear, 'String');
        
        % simulation settings
        Simulation = get(handles.rbSim_Nonlinear, 'Value');
        u          = get(handles.txtU_Sim_Nonlinear, 'String');
        x0         = get(handles.txtX0_Sim_Nonlinear, 'String');
        
        % random simulation settings
        Random_Simulation         = get(handles.rbRSim_Nonlinear, 'Value');
        simOpt_points_random      = get(handles.txtNoP_RSim_Nonlinear, 'String');
        simOpt_fracVert_random    = get(handles.txtFV_RSim_Nonlinear, 'String');
        simOpt_fracInpVert_random = get(handles.txtFIV_RSim_Nonlinear, 'String');
        simOpt_nrConstInp_random  = get(handles.txtNCI_RSim_Nonlinear, 'String');
        
        % RRT simulation settings
        Simulate_RRT          = get(handles.rbSimRRT_Nonlinear, 'Value');
        simOpt_points_RRT     = get(handles.txtNoP_RSim_Nonlinear, 'String');
        Yes                   = get(handles.rbYesEPS_SimRRT_Nonlinear, 'Value');
        No                    = get(handles.rbNoEPS_SimRRT_Nonlinear, 'Value');
        simOpt_stretchFac_RRT = get(handles.txtSF_SimRRT_Nonlinear, 'String');
        
        % plotting settings
        dims = get(handles.listPlots_Nonlinear, 'String');
        intermediate_error_order      = strcmp(get(handles.txtIO_Nonlinear_lin, 'Enable'),'on');
        intermediate_error_order_poly = strcmp(get(handles.txtIO_Nonlinear, 'Enable'),'on');
        
        % reach plot settings
        reach_plot                  = get(handles.checkReachPlot_Nonlinear, 'Value');
        color_contents_reach        = cellstr(get(handles.popColorReach_Nonlinear, 'String'));
        color_reach                 = color_contents_reach{get(handles.popColorReach_Nonlinear, 'Value')};
        edge_color_contents_reach   = cellstr(get(handles.popEdgeColorReach_Nonlinear, 'String'));
        edge_color_reach            = edge_color_contents_reach{get(handles.popEdgeColorReach_Nonlinear, 'Value')};
        
        % initial set plot settings
        initial_plot                = get(handles.checkInitialPlot_Nonlinear, 'Value');
        color_contents_initial      = cellstr(get(handles.popColorInitial_Nonlinear, 'String'));
        color_initial               = color_contents_initial{get(handles.popColorInitial_Nonlinear, 'Value')};
        edge_color_contents_initial = cellstr(get(handles.popEdgeColorInitial_Nonlinear, 'String'));
        edge_color_initial          = edge_color_contents_initial{get(handles.popEdgeColorInitial_Nonlinear, 'Value')};

        % simulation plot settings
        simulation_plot             = get(handles.checkSimulationPlot_Nonlinear, 'Value') && strcmp(get(handles.checkSimulationPlot_Nonlinear, 'Enable'),'on');
        color_contents_simulation   = cellstr(get(handles.popColorSimulation_Nonlinear, 'String'));
        color_simulation            = color_contents_simulation{get(handles.popColorSimulation_Nonlinear, 'Value')};
        line_contents_simulation    = cellstr(get(handles.popLineSimulation_Nonlinear, 'String'));
        line_simulation             = line_contents_simulation{get(handles.popLineSimulation_Nonlinear, 'Value')}; 
        
        % auxiliary function
        Nonlinear_eq = get(handles.rbLoadDynamicEquation_Nonlinear, 'Value') || get(handles.rbEnterDynamicEquation_Nonlinear, 'Value');
        SpaceEx      = get(handles.rbLoadSModel_Nonlinear, 'Value');
        
    case 'Hybrid System'
        
        % name of the function
        func         = 'Hybrid';
        
        % the dimensions and nrOfInputs of the system are needed to make sure 
        % the other parameters are defined correctly 
        dimensions   = handles.dimensions_Hybrid;
        nrOfInputs   = handles.nrOfInputs_Hybrid;
        Sys          = 'hybridSys';
        
        % settings parameters
        tstart       = get(handles.txtTstart_Hybrid, 'String');
        tfinal       = get(handles.txtTfinal_Hybrid, 'String');
        contents     = cellstr(get(handles.popStartLoc_Hybrid, 'String'));
        startLoc     = contents{get(handles.popStartLoc_Hybrid, 'Value')};
        contents     = cellstr(get(handles.popFinalLoc_Hybrid, 'String'));
        finalLoc     = contents{get(handles.popFinalLoc_Hybrid, 'Value')};
        time_step    = get(handles.txtTimeStep_Hybrid, 'String');
        
        % initial and input sets
        R0_handle    = 'R0_Hybrid';
        U_handle     = 'U_Hybrid';
        
        % algorithm settings
        contents     = cellstr(get(handles.popLA_Hybrid, 'String'));
        algorithm    = contents{get(handles.popLA_Hybrid, 'Value')};
        OptError     = get(handles.txtOptionError_Hybrid, 'String');
        zOrder       = get(handles.txtZOrder_Hybrid, 'String');
        taylor_terms = get(handles.txtTaylor_Hybrid, 'String');
        contents     = cellstr(get(handles.popRT_Hybrid, 'String'));
        RT           = contents{get(handles.popRT_Hybrid, 'Value')};
        
        % conservative linearization
        tensor_order_lin        = get(handles.txtConLinTO_Hybrid, 'String');
        intermediate_order_lin  = get(handles.txtIO_Hybrid_lin, 'String');
        error_order_lin         = get(handles.txtEO_Hybrid_lin, 'String');
        
        % conservative polynomilization
        tensor_order_poly       = get(handles.txtConPolTO_Hybrid, 'String');
        intermediate_order_poly = get(handles.txtIO_Hybrid, 'String');
        error_order_poly        = get(handles.txtEO_Hybrid, 'String');
        
        % simulation settings
        Simulation = get(handles.rbSim_Hybrid, 'Value');
        x0         = get(handles.txtX0_Sim_Hybrid, 'String');
        u          = get(handles.txtU_Sim_Hybrid, 'String');
        
        % random simulation settings
        Random_Simulation         = get(handles.rbRSim_Hybrid, 'Value');
        simOpt_points_random      = get(handles.txtNoP_RSim_Hybrid, 'String');
        simOpt_fracVert_random    = get(handles.txtFV_RSim_Hybrid, 'String');
        simOpt_fracInpVert_random = get(handles.txtFIV_RSim_Hybrid, 'String');
        simOpt_nrConstInp_random  = get(handles.txtNCI_RSim_Hybrid, 'String');
        
        % RRT simulation settings
        Simulate_RRT          = get(handles.rbSimRRT_Hybrid, 'Value');
        simOpt_points_RRT     = get(handles.txtNoP_RSim_Hybrid, 'String');
        Yes                   = get(handles.rbYesEPS_SimRRT_Hybrid, 'Value');
        No                    = get(handles.rbNoEPS_SimRRT_Hybrid, 'Value');
        simOpt_stretchFac_RRT = get(handles.txtSF_SimRRT_Hybrid, 'String');
        
        % plotting settings
        dims = get(handles.listPlots_Hybrid, 'String');
        intermediate_error_order      = strcmp(get(handles.txtIO_Hybrid_lin, 'Enable'),'on');
        intermediate_error_order_poly = strcmp(get(handles.txtIO_Hybrid, 'Enable'),'on');
        
        % reach plot settings
        reach_plot                  = get(handles.checkReachPlot_Hybrid, 'Value');
        color_contents_reach        = cellstr(get(handles.popColorReach_Hybrid, 'String'));
        color_reach                 = color_contents_reach{get(handles.popColorReach_Hybrid, 'Value')};
        edge_color_contents_reach   = cellstr(get(handles.popEdgeColorReach_Hybrid, 'String'));
        edge_color_reach            = edge_color_contents_reach{get(handles.popEdgeColorReach_Hybrid, 'Value')};
        
        % initial set plot settings
        initial_plot                = get(handles.checkInitialPlot_Hybrid, 'Value');
        color_contents_initial      = cellstr(get(handles.popColorInitial_Hybrid, 'String'));
        color_initial               = color_contents_initial{get(handles.popColorInitial_Hybrid, 'Value')};
        edge_color_contents_initial = cellstr(get(handles.popEdgeColorInitial_Hybrid, 'String'));
        edge_color_initial          = edge_color_contents_initial{get(handles.popEdgeColorInitial_Hybrid, 'Value')};

        % simulation plot settings
        simulation_plot             = get(handles.checkSimulationPlot_Hybrid, 'Value') && strcmp(get(handles.checkSimulationPlot_Hybrid, 'Enable'),'on');
        color_contents_simulation   = cellstr(get(handles.popColorSimulation_Hybrid, 'String'));
        color_simulation            = color_contents_simulation{get(handles.popColorSimulation_Hybrid, 'Value')};
        line_contents_simulation    = cellstr(get(handles.popLineSimulation_Hybrid, 'String'));
        line_simulation             = line_contents_simulation{get(handles.popLineSimulation_Hybrid, 'Value')};
        
        % for auxiliary function
        SpaceEx = 1; % since Hybrid system has only loading spaceex option
end

% ... Common variables 
% this variable will be used to check if the u in the simulation is
% inside the input set U
u_num = {};

% when uCheck = 0 the default check will be made for both rows and coloum
% otherwise like in the case of adaptive algorithim, only the rows will be
% checked
if strcmp(currentOption, 'Linear System') && strcmp(algorithm, 'adaptive')
    uCheck = 1;
else 
    uCheck = 0;
end

% ... Openeing of the generated file
% This is where the writings on the generated file starts
% To avoid any errors with previously defined variables, the generated file
% is written as a function instead of a script
fprintf(id, '%s, %s\n', '% FILE AUTOMATICALLY GENERATED ON', datetime);

if isempty(func_name)
    fprintf(id, 'function [] = %s()\n\n',func);
else
    fprintf(id, 'function [] = %s()\n\n',erase(func_name,'.m'));
end

fprintf(id, 'disp(''- Running generated file..'')\n\n');

% ... write the first section of the generated file, the system dynamics
fprintf(id, '\n%s\n', '% System Dynamics --------------------------------------------------------- ');

% ... write the system dynamics to the generated file
switch currentOption
    case 'Linear System'
        [err,name] = system_Linear(file, handles,id,err);
    case 'Nonlinear System'
        [err,dimensions,nrOfInputs,name] = system_Nonlinear(file, handles,id,err);
    case 'Hybrid System'
        [err,name] = system_Hybrid(file, handles,id,err);
end
% the err code to prevent the file from running
if err == 1, return, end

% ... write the second section which is the parameter settings to the generated file
fprintf(id, '\n%s\n\n', '% Parameter ---------------------------------------------------------------');

%tstart
if strtrim(tstart) ~= '0'
    fprintf(id, '%s = %s;\n', 'params.tStart', tstart);
end

%tfinal
fprintf(id, '%s = %s;\n', 'params.tFinal', tfinal);

%R0
consPoly = false;
if (strcmp(handles.rbConLin_Hybrid.Visible,'on') && ~get(handles.rbConLin_Hybrid, 'Value')) || ...
   (strcmp(currentOption, 'Nonlinear System') && ~get(handles.rbConLin_Nonlinear, 'Value'))
    consPoly = true;
end
err = R0(handles,id,err,dimensions,R0_handle,consPoly);
if err == 1, return, end

if strcmp(currentOption, 'Hybrid System')
    %startLoc
    fprintf(id, '%s = %s;\n', 'params.startLoc', startLoc);
    
    %finalLoc
    fprintf(id, '%s = %s;\n', 'params.finalLoc', finalLoc);
else
    %u
     [err, u_num] = u_trans(handles,id,err,uTrans_handle,nrOfInputs,tstart,tfinal,time_step,uCheck);
     if err == 1, return, end
end

%U
[err,U_input] = U(handles,id,err,nrOfInputs,U_handle,currentOption);
if err == 1, return, end

% ... write the third section which is the option settings to the generated file
fprintf(id, '\n\n%s\n\n', '% Reachability Analysis ---------------------------------------------------');

fprintf(id, 'disp(''- Running reachability analysis..'')\n\n');

% nonlinear system dynamics or not
nonlin = false;
if strcmp(currentOption, 'Nonlinear System') || ...
   (strcmp(currentOption,'Hybrid System') && ...
    strcmp(handles.rbConLin_Hybrid.Visible,'on'))
    nonlin = true;
end

% linear algorithm    
if ~nonlin
    fprintf(id, '%s = ''%s'';\n', 'options.linAlg', algorithm);
end

% algorithm settings
if nonlin || ~strcmp(algorithm, 'adaptive')
    % time step
    fprintf(id, '%s = %s;\n', 'options.timeStep', time_step);
    % zonotope order
    fprintf(id, '%s = %s;\n', 'options.zonotopeOrder', zOrder);
    % taylor terms
    fprintf(id, '%s = %s;\n', 'options.taylorTerms', taylor_terms);
    % reduction technique
    if ~strcmp(strtrim(RT),'girard')
        fprintf(id, '%s = ''%s'';\n', 'options.reductionTechnique', RT);
    end
end
if strcmp(currentOption, 'Hybrid System') && strcmp(handles.txtGO_Hybrid.Enable,'on')
    % guard order
    gOrder = get(handles.txtGO_Hybrid, 'String');
    fprintf(id, '%s = %s;\n', 'options.guardOrder', gOrder);
end

if ~nonlin && strcmp(algorithm, 'adaptive')
    fprintf(id, '%s = %s;\n', 'options.error', OptError);
end

% conservative linearization or polynomilization
if nonlin
    if (strcmp(handles.rbConLin_Hybrid.Visible,'on') && get(handles.rbConLin_Hybrid, 'Value')) || ...
       (strcmp(currentOption, 'Nonlinear System') && get(handles.rbConLin_Nonlinear, 'Value'))
        % algorithm
        fprintf(id, '%s;\n', 'options.alg = ''lin''');
        % tensor order
        fprintf(id, '%s = %s;\n', 'options.tensorOrder', tensor_order_lin);
        if intermediate_error_order
            % intermidate order
            fprintf(id, '%s = %s;\n', 'options.intermediateOrder', intermediate_order_lin);
            % error order
            fprintf(id, '%s = %s;\n', 'options.errorOrder', error_order_lin);
        end
    elseif get(handles.txtCP_Hybrid, 'Value') || get(handles.rbConPol_Nonlinear , 'Value')
        % algorithm
        fprintf(id, '%s;\n', 'options.alg = ''poly''');
        % tensor order
        fprintf(id, '%s = %s;\n', 'options.tensorOrder', tensor_order_poly);
        if intermediate_error_order_poly
            % intermidate order
            fprintf(id, '%s = %s;\n', 'options.intermediateOrder', intermediate_order_poly);
            % error order
            fprintf(id, '%s = %s;\n', 'options.errorOrder', error_order_poly);
        end
    end
end

% guard intersection: method and enclose (Hybrid system)
if strcmp(currentOption, 'Hybrid System')
    contents = cellstr(get(handles.popMethod_Hybrid, 'String'));
    method = contents{get(handles.popMethod_Hybrid, 'Value')};
    fprintf(id, '%s = ''%s'';\n', 'options.guardIntersect', method);
    % enclose
    if ismember(method,{'polytope','conZonotope','zonoGirard','nondetGuard'})
        if get(handles.checkbox_Hybrid, 'Value')
            if get(handles.checkpca_Hybrid, 'Value')
                if get(handles.checkflow_Hybrid, 'Value')
                    enclose = '{''box'',''pca'',''flow''}';
                else
                    enclose = '{''box'',''pca''}';
                end
            else
                if get(handles.checkflow_Hybrid, 'Value')
                    enclose = '{''box'',''flow''}';
                else
                    enclose = '{''box''}';
                end
            end
        else
            if get(handles.checkpca_Hybrid, 'Value')
                if get(handles.checkflow_Hybrid, 'Value')
                    enclose = '{''pca'',''flow''}';
                else
                    enclose = '{''pca''}';
                end
            else
                if get(handles.checkflow_Hybrid, 'Value')
                    enclose = '{''flow''}';
                else
                    uiwait(msgbox('Unspecified enclose option!', 'Error', 'error', 'modal'))
                    fclose(id);
                    err = 1;
                    return
                end
            end
        end
        fprintf(id, '%s = %s;\n', 'options.enclose', enclose);
    end
end

% reachSet
fprintf(id, '\n%s = %s%s%s;\n\n', 'reachSet', 'reach(',Sys, ' ,params, options)');

% ... write the forth section which is the simulation settings if any of
% the methods were selected
if Simulation
    
    % ... write the option settings to the generated file
    fprintf(id, '\n%s\n\n', '% Simulation --------------------------------------------------------------');
    
    fprintf(id, 'disp(''- Running simulations..'')\n\n');

    % simulation parameters
    if strtrim(tstart) ~= '0'
        fprintf(id, '%s = %s;\n', 'paramsSim.tStart', 'params.tStart');
    end
    fprintf(id, '%s = %s;\n', 'paramsSim.tFinal', 'params.tFinal');
    
    % start and final location for Hybrid system
    if strcmp(currentOption, 'Hybrid System')
        fprintf(id, '%s = %s;\n', 'paramsSim.startLoc', 'params.startLoc');
        fprintf(id, '%s = %s;\n', 'paramsSim.finalLoc', 'params.finalLoc');
    end
    
    % x0 simulation
    err = x0_simulation(x0,id,err,dimensions);
    if err == 1, return, end
    
    % u simulation
    err = u_simulation(u,id,err,nrOfInputs);
    if err == 1, return, end
    
    % check if u is inside U input set
    err = u_check(u,u_num,U_input,id,err);
    if err == 1, return, end
    
    % for Hybrid system the location was also inserted
    if ~strcmp(currentOption, 'Hybrid System')
        fprintf(id, '\n[t,x] = simulate(%s, %s);\n', Sys, 'paramsSim');
        fprintf(id, 'simRes = simResult({x},{t});\n\n');
    else
        fprintf(id, '\n[t,x,loc] = simulate(%s, %s);\n', Sys, 'paramsSim');
        fprintf(id, 'simRes = simResult(x,t,loc);\n\n');
    end
    
elseif Random_Simulation

    % ... write the option settings to the generated file
    fprintf(id, '\n%s\n\n', '% Simulation --------------------------------------------------------------');
    
    fprintf(id, 'disp(''- Running simulations..'')\n\n');
 
    fprintf(id, '%s = %s;\n', 'simOpt.points', simOpt_points_random);
    fprintf(id, '%s = %s;\n', 'simOpt.fracVert', simOpt_fracVert_random);
    fprintf(id, '%s = %s;\n', 'simOpt.fracInpVert', simOpt_fracInpVert_random);
    fprintf(id, '%s = %s;\n', 'simOpt.nrConstInp', simOpt_nrConstInp_random);
    
    fprintf(id, '\n%s = simulateRandom(%s, %s, %s);\n\n', ...
        'simRes', Sys, 'params', 'simOpt');
    
elseif Simulate_RRT
    
    % ... write the option settings to the generated file
    fprintf(id, '\n%s\n\n', '% Simulation --------------------------------------------------------------');
    
    fprintf(id, 'disp(''- Running simulations..'')\n\n');
    
    fprintf(id, '%s = %s;\n', 'simOpt.points', simOpt_points_RRT);
    
    if Yes
        fprintf(id, '%s = %s;\n', 'simOpt.vertSamp', 'boolean(1)');
    elseif No
        fprintf(id, '%s = %s;\n', 'simOpt.vertSamp', 'boolean(0)');
    end
    
    fprintf(id, '%s = %s;\n', 'simOpt.stretchFac', simOpt_stretchFac_RRT);
    fprintf(id, '\n%s = simulateRandom(%s, %s, %s, %s, ''rrt'');\n\n', ...
        'simRes', Sys, 'reachSet', 'params', 'simOpt');
end

% ... write the fifth section which is the plotting settings
fprintf(id, '\n%s\n\n', '% Visualization -----------------------------------------------------------');
    
    fprintf(id, 'disp(''- Plotting..'')\n\n');

gray_color = '[0.5,0.5,0.5]';

if isempty(dims)
    uiwait(msgbox('Missing dimensions to plot', 'Error', 'error', 'modal'))
    fclose(id);
    err = 1;
    return
end

dims = erase(dims,"time, ");

for i = 1:size(dims)
    dims_eva = evalin('base',dims{i});
    if max(dims_eva)> dimensions(1)
        uiwait(msgbox('The dimensions in the plotting list are greater than the system''s dimension! ', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
    end  
end

fprintf(id, '%s\n', '% plot different projections');
fprintf(id, '%s = {', 'dims');
fprintf(id, '%s,', dims{1:size(dims,1)-1});
fprintf(id, '%s', dims{size(dims,1):size(dims,1)});
fprintf(id, '};\n\n%s\n\n', 'for i = 1:length(dims)');
fprintf(id, 'fprintf(''  - Plot #%%i..\\n'',i)\n\n');
fprintf(id, '%s\n', '    figure; hold on; box on');
fprintf(id, '%s = %s\n', '    projDims', 'dims{i};');

if reach_plot
    fprintf(id, '\n%s', '    % plot reachable sets');
    fprintf(id, '\n%s\n', '    if length(projDims) == 1');
    if strcmp(color_reach, 'gray')
        if strcmp(edge_color_reach, 'gray')
             fprintf(id, '        plotOverTime(reachSet, projDims, ''Facecolor'', %s, ''EdgeColor'', %s);\n', gray_color, gray_color);
             fprintf(id, '%s\n', '    else');
             fprintf(id, '        plot(reachSet, projDims,  ''Facecolor'', %s, ''EdgeColor'', %s);\n', gray_color, gray_color);
        else
             fprintf(id, '        plotOverTime(reachSet, projDims,  ''Facecolor'', %s, ''EdgeColor'', ''%s'');\n', gray_color, edge_color_reach);
             fprintf(id, '%s\n', '    else');
             fprintf(id, '        plot(reachSet, projDims,  ''Facecolor'', %s, ''EdgeColor'', ''%s'');\n', gray_color, edge_color_reach);
        end
    else
        if strcmp(edge_color_reach, 'gray')
            fprintf(id, '        plotOverTime(reachSet, projDims,  ''Facecolor'', ''%s'', ''EdgeColor'', %s);\n', color_reach, gray_color);
            fprintf(id, '%s\n', '    else');
            fprintf(id, '        plot(reachSet, projDims,  ''Facecolor'', ''%s'', ''EdgeColor'', %s);\n', color_reach, gray_color);
        else
            fprintf(id, '        plotOverTime(reachSet, projDims,  ''Facecolor'', ''%s'', ''EdgeColor'', ''%s'');\n', color_reach, edge_color_reach);
            fprintf(id, '%s\n', '    else');
            fprintf(id, '        plot(reachSet, projDims,  ''Facecolor'', ''%s'', ''EdgeColor'', ''%s'');\n', color_reach, edge_color_reach);
        end
    end
    fprintf(id, '%s\n', '    end');
end

if initial_plot
    fprintf(id, '\n%s\n', '    % plot initial set');
    fprintf(id, '%s\n', '    if length(projDims) ~= 1');
    if strcmp(color_initial, 'gray')
        if strcmp(edge_color_initial, 'gray')
            fprintf(id, '        plot(params.R0, projDims,  ''Facecolor'', %s, ''EdgeColor'', %s);\n', gray_color, gray_color);
        else
            fprintf(id, '        plot(params.R0, projDims,  ''Facecolor'', %s, ''EdgeColor'', ''%s'');\n', gray_color, edge_color_initial);
        end
    else
        if strcmp(edge_color_initial, 'gray')
            fprintf(id, '        plot(params.R0, projDims,  ''Facecolor'', ''%s'', ''EdgeColor'', %s);\n', color_initial, gray_color);
        else
            fprintf(id, '        plot(params.R0, projDims,  ''Facecolor'', ''%s'', ''EdgeColor'', ''%s'');\n', color_initial, edge_color_initial);
        end
    end
    fprintf(id, '%s\n', '    end');
end

if simulation_plot
    fprintf(id, '\n%s', '    % plot simulation results');
    fprintf(id, '\n%s\n', '    if length(projDims) == 1');
    if strcmp(color_simulation, 'gray')
        fprintf(id, '        plotOverTime(simRes, projDims, ''%s'', ''color'', %s);\n',line_simulation, gray_color); 
        fprintf(id, '%s\n', '    else');
        fprintf(id, '        plot(simRes, projDims, ''%s'', ''color'', %s);\n',line_simulation, gray_color); 
    else
        fprintf(id, '        plotOverTime(simRes, projDims, ''%s'', ''color'', ''%s'');\n',line_simulation, color_simulation); 
        fprintf(id, '%s\n', '    else');
        fprintf(id, '        plot(simRes,projDims, ''%s'', ''color'', ''%s'');\n',line_simulation, color_simulation);  
    end
    fprintf(id, '%s\n', '    end');
end


fprintf(id, '\n%s\n', '    % label plot');
fprintf(id, '%s\n', '    if length(projDims) == 1');
fprintf(id, '        xlabel(''time'')\n' );
fprintf(id, '%s\n', '        ylabel([''x_{'',num2str(projDims(1)),''}''])');
fprintf(id, '%s\n', '    else');
fprintf(id, '%s\n', '        xlabel([''x_{'',num2str(projDims(1)),''}''])');
fprintf(id, '%s\n', '        ylabel([''x_{'',num2str(projDims(2)),''}''])');
fprintf(id, '%s\n\n', '    end');

fprintf(id, 'drawnow %% plot now\n\n');

fprintf(id, '%s\n\n', 'end');

fprintf(id, 'disp(''Done!'')\n\n');

fprintf(id, '%s\n\n', 'end');

% ... write the sixth and final section which is the Auxiliary function
% write SpaceEx function
if SpaceEx
    fprintf(id, '\n%s\n\n', '%% Auxiliary Functions ----------------------------------------------------');
    path = fullfile('models','SpaceExConverted');
    file2 = fullfile(CORAROOT,path,name);
    dinfo = dir(file2);
    dinfo(ismember( {dinfo.name}, {'.', '..'})) = [];
    for i = 1:size(dinfo,1)
        if dinfo(i).isdir
            path_dir = fullfile(file2,dinfo(i).name);
            auxinfo = dir(path_dir);
            auxinfo(ismember( {auxinfo.name}, {'.', '..'})) = [];
            for j = 1:size(auxinfo,1)
                nonlinear_eq = fullfile(path_dir,auxinfo(j).name);
                fID = fopen(nonlinear_eq, 'r');
                fileData = fread(fID);
                fclose(fID);
                fileStr = char(fileData');
                fprintf(id, '%s\n\n', fileStr);
                if ~endsWith(fileStr,'end')
                    fprintf(id, 'end \n\n');
                end
            end
        else
            file = fullfile(file2,dinfo(i).name);
            fID = fopen(file, 'r');
            fileData = fread(fID);
            fclose(fID);
            fileStr = char(fileData');
            fprintf(id, '%s\n', fileStr);
            if ~endsWith(fileStr,'end')
                fprintf(id, 'end \n');
            end
        end
    end
end

% write nonlinear equation
if  strcmp(currentOption, 'Nonlinear System') && Nonlinear_eq
    fprintf(id, '\n%s\n\n', '%% Auxiliary Functions ----------------------------------------------------');
    fID = fopen(handles.dynamicsEq_Nonlinear, 'r');
    fileData = fread(fID);
    fclose(fID);
    fileStr = char(fileData');
    fprintf(id, '%s\n', fileStr);
    if ~endsWith(fileStr,'end')
        fprintf(id, 'end \n\n');
    end
end


fclose(id);
end

% dynamical systems

function [err,name] = system_Linear(~, handles,id,err)
name = {};
if get(handles.rbLoadSModel_Linear, 'Value')
    if isfield(handles, 'spaceExFile')
        name = strsplit(handles.spaceExFile, filesep);
        name = name{end};
        name = strsplit(name, '.');
        name = name{1};
        fprintf(id, '\n%s = %s;\n\n', 'linSys', strcat(name, '()'));
    else
        uiwait(msgbox('Undefined spaceEx File', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
    end
elseif get(handles.rbCreateMat_Linear, 'Value')
    if get(handles.checkA_Linear, 'Value')
        if isfield(handles, 'A')
            mat_str = handles.A;
            A = get(handles.textA_Linear, 'String');
            if ~isempty(mat_str)
                fprintf(id, '\n%s = %s;\n', 'A', A);
                mats_str = {'A'};
            else
                uiwait(msgbox('Matrix A must be defined', 'Error', 'error', 'modal'));
                fclose(id);
                err = 1;
                return
            end
        else
            uiwait(msgbox('Matrix A must be defined', 'Error', 'error', 'modal'));
            fclose(id);
            err = 1;
            return
        end
    else
        uiwait(msgbox('Matrix A must be defined', 'Error', 'error', 'modal'));
        fclose(id);
        err = 1;
        return
    end
    
    if get(handles.checkB_Linear, 'Value')
        if isfield(handles, 'B')
            mat_str = handles.B;
            B = get(handles.textB_Linear, 'String');
            temp = evalin('base',B);
            A = evalin('base',A);
            if isscalar(temp)
                handles.nrOfInputs_Linear = size(A,1);
            else
                handles.nrOfInputs_Linear = size(temp,2);
            end
            if ~isempty(mat_str)
                fprintf(id, '%s = %s;\n', 'B', B);
                mats_str = join({mats_str{1},'B'}, ',');
            else
                uiwait(msgbox('Matrix B must be defined', 'Error', 'error', 'modal'))
                fclose(id);
                err = 1;
                return
            end
        else
            uiwait(msgbox('Matrix B must be defined', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
    else
        uiwait(msgbox('Matrix B must be defined', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
    end
    
    if get(handles.checkC_Linear, 'Value')
        if isfield(handles, 'C')
            mat_str = handles.C;
            c = get(handles.textC_Linear, 'String');
            if ~isempty(mat_str)
                fprintf(id, '%s = %s;\n', 'c', c);
                mats_str = join({mats_str{1},'c'}, ',');
            else
                uiwait(msgbox('Matrix c is not defined', 'Error', 'error', 'modal'))
                fclose(id);
                err = 1;
                return
            end
        else
            uiwait(msgbox('Matrix c is not defined', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
    end
    
    if get(handles.checkD_Linear, 'Value')
        if isfield(handles, 'D')
            mat_str = handles.D;
            C = get(handles.textD_Linear, 'String');
            if ~isempty(mat_str)
                fprintf(id, '%s = %s;\n', 'C', C);
                mats_str = join({mats_str{1},'C'}, ',');
            else
                uiwait(msgbox('Matrix C is not defined', 'Error', 'error', 'modal'))
                fclose(id);
                err = 1;
                return
            end
        else
            uiwait(msgbox('Matrix C is not defined', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
    end
    
    if get(handles.checkE_Linear, 'Value')
        if isfield(handles, 'E')
            mat_str = handles.E;
            D = get(handles.textE_Linear, 'String');
            if ~isempty(mat_str)
                fprintf(id, '%s = %s;\n', 'D', D);
                mats_str = join({mats_str{1},'D'}, ',');
            else
                uiwait(msgbox('Matrix D is not defined', 'Error', 'error', 'modal'))
                fclose(id);
                err = 1;
                return
            end
        else
            uiwait(msgbox('Matrix D is not defined', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
    end
    
    if get(handles.checkF_Linear, 'Value')
        if isfield(handles, 'F')
            mat_str = handles.F;
            K = get(handles.textF_Linear, 'String');
            if ~isempty(mat_str)
                fprintf(id, '%s = %s;\n', 'K', K);
                mats_str = join({mats_str{1},'K'}, ',');
            else
                uiwait(msgbox('Matrix K is not defined', 'Error', 'error', 'modal'))
                fclose(id);
                err = 1;
                return
            end
        else
            uiwait(msgbox('Matrix K is not defined', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
    end
    
    fprintf(id, '\n%s = linearSys(''%s'', %s);\n\n', 'linSys', 'linearSys', ...
        mats_str{1});
else
    uiwait(msgbox('Linear system is not defined. Please either load a SpaceEx file or create the system using matrices', ...
        'Error', 'error', 'modal'))
    fclose(id);
    err = 1;
    return
end
end

function [err,dimensions,nrOfInputs,name] = system_Nonlinear(~, handles,id,err)
dimensions = {};
nrOfInputs = {};
name = {};
if get(handles.rbLoadSModel_Nonlinear, 'Value')
    handles.dynamicsEq_Nonlinear = [];
    if isfield(handles, 'spaceExFile_Nonlinear')
        name = strsplit(handles.spaceExFile_Nonlinear, filesep);
        name = name{end};
        name = strsplit(name, '.');
        name = name{1};
        dimensions = handles.dimensions_Nonlinear;
        nrOfInputs = handles.nrOfInputs_Nonlinear;
        
        fprintf(id, '\nnonlinSys = %s;\n\n', strcat(name, '()'));
    else
        uiwait(msgbox('Undefined spaceEx File', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
    end
elseif get(handles.rbLoadDynamicEquation_Nonlinear, 'Value') || get(handles.rbEnterDynamicEquation_Nonlinear, 'Value')
    handles.spaceExFile_Nonlinear = [];
    if isfield(handles, 'dynamicsEq_Nonlinear') && handles.dynamicsEq_Nonlinear_error
        dynamics_equation = handles.dynamicsEq_Nonlinear;
        dynamics_equation = strsplit(dynamics_equation, filesep);
        dynamics_equation = dynamics_equation{end};
        dynamics_equation = strsplit(dynamics_equation,'.');
        dynamics_equation = dynamics_equation{1};
        dynamicsEq_handle_str = ['@', dynamics_equation];
        dynamicsEq_handle = eval('base', dynamicsEq_handle_str);
        
        %extract num_dimensions, num_inputs
        [temp,dimensions]  = inputArgsLength_app(dynamicsEq_handle,2);
        nrOfInputs = max(1,temp(2));
        
        fprintf(id, '\nnonlinSys = nonlinearSys(%s);\n\n',dynamicsEq_handle_str);
    elseif ~handles.dynamicsEq_Nonlinear_error
        uiwait(msgbox('Error in the dynamical equation', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
    else
        uiwait(msgbox('Undefined dynamic equation', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
    end
else
    uiwait(msgbox('Nonlinear system is not defined. Please either load a SpaceEx file, a dynamic equation file or enter dynamics equation', ...
        'Error', 'error', 'modal'))
    fclose(id);
    err = 1;
    return
end

end

function [err,name] = system_Hybrid(~, handles,id,err)
name = {};
if isfield(handles, 'spaceExFile_Hybrid')
    name = strsplit(handles.spaceExFile_Hybrid, filesep);
    name = name{end};
    name = strsplit(name, '.');
    name = name{1};
    fprintf(id, '\nhybridSys = %s;\n\n', strcat(name, '()'));
else
    uiwait(msgbox('Undefined spaceEx File', 'Error', 'error', 'modal'))
    fclose(id);
    err = 1;
    return
end

end

% initial and input sets

function [err] = R0(handles,id,err,dimensions,R0_handle,consPoly)

if isfield(handles, R0_handle)
    R0 = eval('base',['handles.',R0_handle]);
    if ~isempty(R0)
        type = R0.type;
        center_str = R0.center;
        width_str = R0.width;
        LB = evalin('base', center_str);
        UB = evalin('base', width_str);
        R0_size = size(LB,1);
        if dimensions ~= R0_size
            uiwait(msgbox('dimension of R0 is not the same as the system', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
        if strcmp(type, 'interval')
            [mrows, ncols] = size(UB);
            if mrows == 1 && ncols == 1
                if ~consPoly
                    fprintf(id, '%s = zonotope(interval(%d, %d));\n', 'params.R0',LB,UB);
                else
                    fprintf(id, '%s = polyZonotope(zonotope(interval(%d, %d)));\n', 'params.R0',LB,UB);
                end
            else
                outputstr = ['%.' num2str(mrows) 'f '];
                outputstr = repmat(outputstr, 1, ncols);
                if ~consPoly
                    fprintf(id,'%s = zonotope(interval([', 'params.R0');
                else
                    fprintf(id,'%s = polyZonotope(zonotope(interval([', 'params.R0');
                end
                for i = 1:mrows
                    fprintf(id, outputstr, LB(i,:).');
                    if i ~= mrows
                        fprintf(id, ';');
                    end
                end
                fprintf(id, '],[');
                for i = 1:mrows
                    fprintf(id, outputstr, UB(i,:).');
                    if i ~= mrows
                        fprintf(id, ';');
                    end
                end
                if ~consPoly
                    fprintf(id, ']));\n');
                else
                    fprintf(id, '])));\n');
                end
            end
        else
            if ~consPoly
                fprintf(id, '%s = zonotope(%s, %s);\n', 'params.R0', center_str, width_str);
            else
                fprintf(id, '%s = polyZonotope(zonotope(%s, %s));\n', 'params.R0', center_str, width_str);
            end
        end
    else
        uiwait(msgbox('Variable R0 in settings must be defined', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
    end
else
    uiwait(msgbox('Variable R0 in settings must be defined', 'Error', 'error', 'modal'))
    fclose(id);
    err = 1;
    return
end
end

function [err,U_input] = U(handles,id,err,nrOfInputs,U_handle,currentOption)
% the U_input was defined here with as an empty array to prevent having an 
% in the output argument in case U was empty
U_input = [];
if isfield(handles, U_handle)
    U = eval('base',['handles.',U_handle]);
    if ~isempty(U)
        if strcmp(currentOption, 'Hybrid System')
            if ~isempty(U{1,1})
                type = U{1,1}.type;
                center_str = U{1,1}.center;
                width_str = U{1,1}.width;
                p =0;
                [err,U_input] = U_function(id,center_str,width_str,p,nrOfInputs,type,err);
                if err == 1, return, end
            else
                for p = 2:size(U,2)
                    type = U{1,p}.type;
                    center_str = U{1,p}.center;
                    width_str = U{1,p}.width;
                    [err,U_input] = U_function(id,center_str,width_str,p,nrOfInputs,type,err);
                    if err == 1, return, end
                end
                
            end
        else
            type = U.type;
            center_str = U.center;
            width_str = U.width;
            p =0;
            [err,U_input] = U_function(id,center_str,width_str,p,nrOfInputs,type,err);
            if err == 1, return, end
        end
    else
        uiwait(msgbox('Variable U in settings must be defined', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
    end
else
    uiwait(msgbox('Variable U in settings must be defined', 'Error', 'error', 'modal'))
    fclose(id);
    err = 1;
    return
end
end

function [err,U_input] = U_function(id,center_str,width_str,p,nrOfInputs,type,err)
% the U_input was defined here with as an empty array to prevent having an 
% in the output argument in case U was empty
U_input = [];
LB = evalin('base', center_str);
UB = evalin('base', width_str);
U_size = size(LB,1);
if nrOfInputs ~= U_size
    uiwait(msgbox('dimension of U is not the same as the input state of the system', 'Error', 'error', 'modal'))
    fclose(id);
    err = 1;
    return
end
if strcmp(type, 'interval')
    [mrows, ncols] = size(UB);
    if mrows == 1 && ncols == 1
        if p > 1
            fprintf(id, '%s{%d} = zonotope(interval(%d, %d));\n', 'params.U',p-1,LB,UB);
        else
            fprintf(id, '%s = zonotope(interval(%d, %d));\n', 'params.U',LB,UB);
        end  
    else
        outputstr = ['%.' num2str(mrows) 'f '];
        outputstr = repmat(outputstr, 1, ncols);
        fprintf(id,'%s = zonotope(interval([', 'params.U');
        for i = 1:mrows
            fprintf(id, outputstr, LB(i,:).');
            if i ~= mrows
                fprintf(id, ';');
            end
        end
        fprintf(id, '],[');
        for i = 1:mrows
            fprintf(id, outputstr, UB(i,:).');
            if i ~= mrows
                fprintf(id, ';');
            end
        end
        fprintf(id, ']));\n');
    end
    if p > 1
        U_input{p} = interval(LB,UB);
    else
        U_input = interval(LB,UB);
    end
else
    if p > 1
        U_input{p} = zonotope(LB,UB);
        fprintf(id, '%s{%d} = zonotope(%s, %s);\n', 'params.U',p-1, center_str, width_str);
    else
        U_input = zonotope(LB,UB);
        fprintf(id, '%s = zonotope(%s, %s);\n', 'params.U', center_str, width_str);
    end   
end
end

function [err, u_num] = u_trans(handles,id,err,uTrans_handle,nrOfInputs,tstart,tfinal,time_step,uCheck)
% the u_num was defined here with as an empty array to prevent having an 
% in the output argument in case U was empty
u_num = [];
if isfield(handles,uTrans_handle )
    u = ['handles.',uTrans_handle];
    u_trans = eval('base',u);
    if ~isempty(u_trans)
        u_num = evalin('base', u_trans);
        if nrOfInputs ~= size(u_num,1)
            uiwait(msgbox('Number of rows of u is not the same as the input state of the system', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
        if ~uCheck
            tstart_num = evalin('base', tstart);
            tfinal_num = evalin('base', tfinal);
            time_step_num  = evalin('base', time_step);
            if (tfinal_num-tstart_num)/time_step_num ~= size(u_num,2)
                uiwait(msgbox('Number of columns of u is not the same as the (final time - start time)/time step', 'Error', 'error', 'modal'))
                fclose(id);
                err = 1;
                return
            end
        end
        fprintf(id, '%s = %s;\n', 'params.u', u_trans);
    end
end
end

% simulation settings

function [err] = x0_simulation(x0,id,err,dimensions)
if ~isempty(x0)
        x0_size = evalin('base', x0);
        if dimensions ~= length(x0_size)
            uiwait(msgbox('dimension of x0 is not the same as the system', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
        fprintf(id, '%s = %s;\n', 'paramsSim.x0', x0);
    else
        uiwait(msgbox('Variable x0 in simulation is not defined', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
end
end

function [err] = u_simulation(u,id,err,nrOfInputs)
if ~isempty(u)
        u_size = evalin('base', u);
        if nrOfInputs ~= size(u_size,1)
            uiwait(msgbox('Number of rows of u in the simulation is not the same as the input state of the system', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
        fprintf(id, '%s = %s;\n', 'paramsSim.u', u);
    else
        uiwait(msgbox('Variable u in simulation is not defined', 'Error', 'error', 'modal'))
        fclose(id);
        err = 1;
        return
end
end

function [err] = u_check(u,u_num,U_input,id,err)
u_size = evalin('base', u);
k = 0 ;
if isempty(u_num)
        for i = 1:size(u_size,2)
            if size(U_input,2) == 1
                if ~contains(U_input,u_size(:,i))
                    uiwait(msgbox('The value of u in simulation is not inside the input set (U)', 'Error', 'error', 'modal'))
                    fclose(id);
                    err = 1;
                    return
                end
            else
                for j = 1:size(U_input,2)
                    if ~isempty(U_input{j})
                        if ~contains(U_input{j},u_size(:,i))
                            k = k + 1;
                        end
                    end
                end
                if k == size(U_input,2)-1
                    uiwait(msgbox('The value of u in simulation is not inside the input set (U)', 'Error', 'error', 'modal'))
                    fclose(id);
                    err = 1;
                    return
                end
            end
        end
elseif size(u_num,2) == size(u_size,2)
    for i = 1:size(u_size,2)
        if ~contains(U_input+u_num(:,i),u_size(:,i))
            uiwait(msgbox('The value of u in simulation is not inside the input set (U)', 'Error', 'error', 'modal'))
            fclose(id);
            err = 1;
            return
        end
    end  
else
    time = 0:1/size(u_num,2):1;
    timeSim = 0:1/size(u_size,2):1;
    for i = 1:size(u_num,2)
        for j = 1:size(u_size,2)
            inf1 = time(i); sup1 = time(i+1); inf2 = timeSim(j); sup2 = timeSim(j+1);
            if ((sup1 > sup2 + eps) && (inf1 + eps < sup2)) || ((inf1 + eps < inf2) && (sup1 > inf2 + eps))
                if ~contains(U_input+u_num(:,i),u_size(:,j))
                    uiwait(msgbox('The value of u in simulation is not inside the input set (U)', 'Error', 'error', 'modal'))
                    fclose(id);
                    err = 1;
                    return
                end
            end
        end
    end
end
end
