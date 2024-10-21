function [R, eval] = validateReach(obj, configs, check_contain, plot_settings)
% validateReach - evaluate the given test case, i.e., compute the reachable
%   set of the systems in configs, check containment of all measurements,
%   and plot the results
%
% Syntax:
%    R = validateReach(obj, configs)
%    R = validateReach(obj, configs, check_contain)
%    R = validateReach(obj, configs, check_contain, plot_settings)
%
% Inputs:
%    obj - testCase object
%    configs - cell array containing the configurations, i.e.,
%        {i}.sys - dynamical system
%        {i}.params - reachability parameters
%        {i}.options - reachability options
%        {i}.name [optional] - dynamical system
%    check_contain - boolean to check containment of the measurements (only if defined)
%    plot_settings - settings for generating plots (only if defined)
%        .plot_Yp - boolean for plotting the partial reachable set Y_p
%           default: false
%        .k_plot - array containing the time steps which will be plotted
%           default: 1:n_k
%        .dims - n x 2 array containing the output dimensions which will be plotted
%           default: [1 2]
%        .path_fig - path for saving the figure (only saved if specified)
%           default: "" (not saved)
%        .k_save - array containing the time steps which will be saved
%           default: 1:n_k
%        .name - title for the plot
%           default: "Reachability for y_1 and y_2"
%        .s_val - start index for the samples from the validation testcases
%           default: Inf
%        .colors - array of the line colors
%           default: [gray; colors from colormap('lines')]
%        .lines - array of the line types
%           default: ["-", "-", "--", "-.", ":"]
%        .linewidths - array of the line widths
%           default: ["-", "-", "--", "-.", ":"]
%
% Outputs:
%    R - reachable sets
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: testCase

% Authors:       Laura Luetzow
% Written:       01-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

R = cell(length(configs),1);
R_p = cell(length(configs),1);

% Evaluate each system in configs -----------------------------------------

eval.num_out = zeros(length(configs),1);
eval.num_in = zeros(length(configs),1);
eval.size_R = zeros(length(configs),1);
for i = 1: length(configs)
    % 1. compute the reachable set R

    if ~isfield(configs{i}.params,'U')
        for i_k = 1: size(obj.u,1)
            R{i}{i_k} = [];
        end
        continue
    end
    y = obj.y;
    params.U = configs{i}.params.U;

    if size(obj.u,3) > 1 
        % plot y_a instead of y
        configs{i}.obj = obj.compute_ya(configs{i}.sys);
        y = configs{i}.obj.y_a;
        params.R0 = configs{i}.params.R0;
        params.u = zeros(size(obj.u,2),size(obj.u,1));
    elseif isa(configs{i}.sys, 'nonlinearARX') || isa(configs{i}.sys, 'linearARX')
        % input-output model
        obj2 = setInitialStateToMeas(obj,configs{i}.sys.n_p, 1e-2);
        if length(obj2) > 1
            throw(CORAerror('CORA:specialError','Testcase is not valid'))
            % measurement trajectories must be identical for time steps k=1,..,p
        end
        if sum(configs{i}.params.R0.G, 'all') > 1e-6
            throw(CORAerror('CORA:notSupported','Uncertian R0 not implemented.'))
        end
        if isa(configs{i}.sys, 'nonlinearARX')
            params.R0 = configs{i}.params.R0 + obj2{1}.initialState;
        else
            params.y0 = reshape(center(configs{i}.params.R0) + obj2{1}.initialState, configs{i}.sys.nrOfOutputs, []);
            params.R0 = configs{i}.params.R0 + obj.initialState;
        end
        params.u = obj.u';
    else
        params.u = obj.u';
        params.R0 = configs{i}.params.R0 + obj.initialState;
    end

    diff_u = configs{i}.sys.nrOfInputs - size(params.u,1);
    if diff_u ~= 0
        params.u = cat(1, params.u, zeros(diff_u, size(params.u,2), size(params.u,3)));
    end

    % compute reachable set
    if isfield(configs{i}.options, 'cs')
        options = rmfield(configs{i}.options,'cs');
    else
        options = configs{i}.options;
    end    
    n_k = size(params.u,2);
    params.tFinal = configs{i}.sys.dt * (n_k-1);
    R{i} = reach(configs{i}.sys, params, options);
    R{i} = R{i}.timePoint.set;

    % 2. check containment
    if nargin > 2 && check_contain
        for k=1:n_k
            for s=1:size(y,3)
                if ~isa(R{i}{k},'contSet') || representsa(R{i}{k},'emptySet') || any(isnan(y(k,:,s)))
                    % not valid
                elseif (isa(R{i}{k},'zonotope') && (isempty(R{i}{k}.generators) || ...
                        sum(abs((R{i}{k}.generators)),'all') == 0) && ...
                        sum(abs(R{i}{k}.center - y(k,:,s)'),'all') >1e-6) || ...
                        ~contains(R{i}{k}, y(k,:,s)', 'approx', 1e-6)
                    eval.num_out(i) = eval.num_out(i) + 1;
                else
                    eval.num_in(i) = eval.num_in(i) + 1;
                end
            end
            eval.size_R(i) = eval.size_R(i) + sum(abs((R{i}{k}.generators)),'all');
        end
    end
    if isa(configs{1}.sys, 'nonlinearSysDT') && isa(configs{i}.sys, 'nonlinearARX')
        p = configs{1}.sys.nrOfStates;
        R{i} = R{i}(p+1:end);
        y = obj.y(p+1:end,:,:);
    end

    % compute the partial reachable set R_p due to linearization errors
    if nargin >= 4 && isfield(plot_settings,'plot_Yp') && plot_settings.plot_Yp
        p_GO = computeGO(configs{i}.sys, params.R0.c, params.U.c + obj.u', n_k);
        Y_p = cell(n_k,1);
        for k=1:n_k
            Y_p{k} = p_GO.y(:,k) + p_GO.C{k} * ...
                (configs{i}.params.R0-center(configs{i}.params.R0));
            for j = 1:k
                Y_p{k} = Y_p{k} + p_GO.D{k,j} * ...
                    (configs{i}.params.U-center(configs{i}.params.U));
            end
        end
        R_p{i} = Y_p;
    end
end

if nargin < 4 || (isfield(plot_settings,'k_plot') && isempty(plot_settings.k_plot))
    return
end

% Plot the reachable sets and the identification data ---------------------

% plotting specifications
if ~isfield(plot_settings, 'lines')
    line = repmat(["-", "-", "--", "-.", ":"], 1, ceil(length(configs)/4));
else
    line = plot_settings.lines;
end
if ~isfield(plot_settings, 'colors')
    c = colormap("lines");
else
    c = plot_settings.colors;
end
if ~isfield(plot_settings, 'linewidths')
    lw = 0.1 * ones(length(configs));
else
    lw = plot_settings.linewidths;
end
if ~isfield(plot_settings,'k_plot')
    plot_settings.k_plot = 1:n_k;
end
if ~isfield(plot_settings,'dims')
    plot_settings.dims = [1 2];
end
num_row = ceil(sqrt(length(plot_settings.k_plot)));
num_col = ceil(length(plot_settings.k_plot)/num_row);

% create a figure for each dimension pair specified in plot_settings.dims
for i_dim = 1 : size(plot_settings.dims,1)
    figure('units','normalized','outerposition',[0 0 1 1]);
    dim = plot_settings.dims(i_dim,:);
    if ~isfield(plot_settings,'name')
        plot_settings.name = sprintf('Reachability for y_%d and y_%d', dim(1), dim(2));
    end

    % subplot for each time step specified in plot_settings.k_plot
    for i_k = 1:length(plot_settings.k_plot)
        subplot(num_col, num_row, i_k); hold on;
        k = plot_settings.k_plot(i_k);

        % plot the reachable set for each system in configs
        for i=1:length(configs)
            if isfield(configs{i},'name') && ~isempty(configs{i}.name)
                name = configs{i}.name;
            else
                name = sprintf('config%d',i);
            end
            if contains(name, "true")
                col = [0.9 0.9 0.9];
                alpha = 0.8;
            else
                col = c(i,:);
                alpha = 0;
            end
            if size(obj.u,3) > 1
                name_Y = "$\mathcal{Y}_{a,\displaystyle \mathrm{" ...
                    + name + "}}$";
            else
                name_Y = "$\mathcal{Y}_{\displaystyle \mathrm{" ...
                    + name + "}}$";
            end
            if ~isa(R{i}{k}, 'contSet')
                R_ik = zonotope(R{i}{k});
            else 
                R_ik = R{i}{k};
            end
            if isempty(R_ik.generators) || ...
                    sum(abs((R_ik.generators)),'all') == 0
                if isempty(R_ik)
                    continue
                end
                plot(R_ik,[dim(1),dim(2)], 'FaceColor', col, ...
                    'Marker', 'o', 'LineWidth', lw(i), 'DisplayName', name_Y);
            else
                if sum(abs(R_ik.generators),'all') < 10^8 || ~contains(name_Y,'gp')
                    plot(R_ik,[dim(1),dim(2)], 'FaceColor', col, ...
                        'FaceAlpha',alpha, 'LineWidth', lw(i), ...
                        'LineStyle', line(i), 'DisplayName', name_Y);
                end
            end

            % maybe plot the partial reachable set Y_p
            if ~isempty(R_p{i})
                name_Yp = "$\mathcal{Y}_{\displaystyle p,\mathrm{" ...
                    + name + "}}$";
                if isempty(R_p{i}{k}.generators) || ...
                        sum(abs((R_p{i}{k}.generators)),'all') == 0
                    plot(R_p{i}{k},[dim(1),dim(2)], 'FaceColor', col, ...
                        'Marker', 'o', 'LineWidth', 0.5, 'DisplayName', name_Yp);
                else
                    plot(R_p{i}{k},[dim(1),dim(2)],'FaceColor', col, ...
                        'FaceAlpha',alpha, 'LineWidth', 0.5, ...
                        'LineStyle', line(i), 'DisplayName', name_Yp);
                end
            end
            if size(obj.u,3) > 1
                % plot for each system the computed y_a
                plot(squeeze(configs{i}.obj.y_a(k,dim(1),:)), squeeze(configs{i}.obj.y_a(k,dim(2),:)), ...
                    'x', 'Color', col, 'HandleVisibility', 'off');
            end
        end

        % plot the identification data
        if size(obj.u,3) <= 1
            y1 = squeeze(y(k,dim(1),:));
            y2 = squeeze(y(k,dim(2),:));
            name_y = '$y^{\displaystyle(s)}$';
        else
            y1 = squeeze(configs{i}.obj.y_a(k,dim(1),:));
            y2 = squeeze(configs{i}.obj.y_a(k,dim(2),:));
            name_y = "$y_{\displaystyle a}^{\displaystyle(s)}$";
        end
        
        if isfield(plot_settings,'s_val') 
            % plot validation test cases in gray
            name_yval = "$y_{\displaystyle a,val}^{\displaystyle(s)}$";    
            plot(y1(plot_settings.s_val:end), y2(plot_settings.s_val:end), ...
                'x', 'Color', [0.8 0.8 0.8], 'DisplayName', name_yval);
            plot(y1(1:plot_settings.s_val-1), y2(1:plot_settings.s_val-1), 'kx', 'DisplayName', name_y);
        else
            plot(y1, y2, 'kx', 'DisplayName', name_y);
        end

        % create legend
        if k==plot_settings.k_plot(1)
            lg = legend('location', 'best', 'AutoUpdate', ...
                'off', 'Interpreter','latex', 'FontSize',10);
            set(lg.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData',uint8(255*[1;1;1;.9]));
        end

        % adapt format and save specified plots
        fontname(gca,"times");
        if k==1
            text_k = "k=n_p" + "\newline\newline  " ;
        else
            text_k = sprintf("k=n_p+%d", k-1) + "\newline\newline    " ;
        end
        xlabel(sprintf('y_%d', dim(1)));
        if isa(configs{1}.sys, 'linearSysDT')
            ylabel(text_k + sprintf('y_%d', dim(2)));
        else
            ylabel(sprintf('y_%d', dim(2)));
        end
        if ~isfield(plot_settings,'path_fig')
            title("Reachable Output Sets at Timestep " + k);
        elseif ~isfield(plot_settings,'k_save') || any(plot_settings.k_save==k)
            legend('off')
            xl = xlim;
            yl = ylim;
            xlim([xl(1)-0.05*(xl(2)-xl(1)), xl(2)+0.05*(xl(2)-xl(1))])
            ylim([yl(1)-0.05*(yl(2)-yl(1)), yl(2)+0.05*(yl(2)-yl(1))])
            exportgraphics(gca,plot_settings.path_fig+"_k"+k+".pdf",...
             'ContentType','vector')
        end
    end
    sgtitle(plot_settings.name)
end
end

% ------------------------------ END OF CODE ------------------------------
