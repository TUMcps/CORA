function [R, eval] = validateReach(traj, configs, check_contain, plot_settings)
% validateReach - evaluate the given trajectory, i.e., compute the reachable
%   set of the systems in configs, check containment of all measurements,
%   and plot the results
%
% Syntax:
%    R = validateReach(obj, configs)
%    R = validateReach(obj, configs, check_contain)
%    R = validateReach(obj, configs, check_contain, plot_settings)
%
% Inputs:
%    traj - trajectory object
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
%        .s_val - start index for the samples from the validation trajectories
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
% See also: trajectory

% Authors:       Laura Luetzow
% Written:       28-August-2025
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
        for i_k = 1: size(traj.u,1)
            R{i}{i_k} = [];
        end
        continue
    end
    [params, y, plot_ya] = aux_setParams(configs{i}, traj);
    configs{i}.y = y;

    % call reach
    if isfield(configs{i}.options, 'cs')
        options = rmfield(configs{i}.options,'cs');
    else
        options = configs{i}.options;
    end    
    R{i} = reach(configs{i}.sys, params, options);
    R{i} = R{i}.timePoint.set;

    % 2. check containment
    if nargin > 2 && check_contain
        [eval.num_out(i), eval.num_in(i), eval.size_R(i)] = aux_checkContain(R{i}, y);
    end
    % if isa(configs{1}.sys, 'nonlinearSysDT') && isa(configs{i}.sys, 'nonlinearARX')
    %     % nonlinear model
    %     p = configs{1}.sys.nrOfDims;
    %     R{i} = R{i}(p+1:end);
    %     y = traj.y(p+1:end,:,:);
    % end

    % compute the partial reachable set R_p due to linearization errors
    if nargin >= 4 && isfield(plot_settings,'plot_Yp') && plot_settings.plot_Yp
        R_p{i} = aux_computeRp(configs{i},params,traj);
    end
end

if nargin < 4 || (isfield(plot_settings,'k_plot') && isempty(plot_settings.k_plot))
    % don't create plots
    return
end

% Plot the reachable sets and the identification data ---------------------
aux_plot(configs,plot_settings,y,plot_ya,R,R_p)

end


% Auxiliary functions -----------------------------------------------------

function [params, y, plot_ya] = aux_setParams(config_i,traj)
% set parameters for reachability analysis


% set input and initial state set
if length(traj) > 1 && (isa(config_i.sys, 'linearSysDT') ...
        || isa(config_i.sys, 'linearARX'))
    % linear models: visualize all trajectories in one plot using
    % superposition principle (plot y_a instead of y)
    plot_ya = true;
    y = traj.computeOutputDev(config_i.sys);
    params.R0 = config_i.params.R0;
    params.u = zeros(size(traj(1).u,1),size(y,2));
elseif length(traj) > 1
    throw(CORAerror('CORA:specialError',['For nonlinear systems, call' ...
        'validateReach with a single trajectory object traj(i).']));
else
    plot_ya = false;
    y = traj.y;
    params.u = traj.u;
    params.R0 = config_i.params.R0 + traj.x(:,1,1);
end
n_k = size(y,2);
params.tFinal = config_i.sys.dt * (n_k-1);
params.U = config_i.params.U;
end

function [num_out, num_in, size_R] = aux_checkContain(R_i, y)
% check containment of the measurements in the reachablse sets
num_out = 0;
num_in = 0;
size_R = 0;
n_k = size(y,2);
for k=1:n_k
    for s=1:size(y,3)
        if ~isa(R_i{k},'contSet') || representsa(R_i{k},'emptySet') || any(isnan(y(:,k,s)))
            % not valid
        elseif (isa(R_i{k},'zonotope') && (isempty(R_i{k}.generators) || ...
                sum(abs((R_i{k}.generators)),'all') == 0) && ...
                sum(abs(R_i{k}.center - y(:,k,s)),'all') >1e-6) || ...
                ~contains(R_i{k}, y(:,k,s), 'approx', 1e-6)
            num_out = num_out + 1;
        else
            num_in = num_in + 1;
        end
    end
    size_R = size_R + sum(abs((R_i{k}.generators)),'all');
end
end

function R_p_i = aux_computeRp(config_i,params,traj)
% compute the partial reachable set R_p due to linearization errors
n_k = size(params.u,2);
p_GO = computeGO(config_i.sys, params.R0.c, params.U.c + traj.u', n_k);
R_p_i = cell(n_k,1);
for k=1:n_k
    R_p_i{k} = p_GO.y(:,k) + p_GO.C{k} * ...
        (config_i.params.R0-center(config_i.params.R0));
    for j = 1:k
        R_p_i{k} = R_p_i{k} + p_GO.D{k,j} * ...
            (config_i.params.U-center(config_i.params.U));
    end
end
end

function aux_plot(configs,plot_settings,y,plot_ya,R,R_p)
% plotting specifications
if ~isfield(plot_settings,'k_plot')
    plot_settings.k_plot = 1:size(y,2);
end

if ~isfield(plot_settings,'dims')
    plot_settings.dims = [1 2];
end

% set line types
if ~isfield(plot_settings, 'lines')
    line = repmat(["-", "-", "--", "-.", ":"], 1, ceil(length(configs)/4));
else
    line = plot_settings.lines;
end

% set line colors
if ~isfield(plot_settings, 'colors')
    c = colormap("lines");
else
    c = plot_settings.colors;
end

% set line width
if ~isfield(plot_settings, 'linewidths')
    lw = 0.1 * ones(length(configs));
else
    lw = plot_settings.linewidths;
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
            % get name
            if isfield(configs{i},'name') && ~isempty(configs{i}.name)
                name = configs{i}.name;
            else
                name = sprintf('config%d',i);
            end

            % get display name
            name_Y = "$\mathcal{Y}_{\displaystyle \mathrm{" ...
                + name + "}}$";
            name_y = '$y^{\displaystyle(s)}$';
            name_yval = "$y_{\displaystyle val}^{\displaystyle(s)}$";

            % get color and alpha
            if contains(name, "true")
                col = [0.9 0.9 0.9];
                alpha = 0.8;
            else
                col = c(i,:);
                alpha = 0;
            end
            % build as set
            if length(R{i}) >= k
                R_ik = R{i}{k};
                if ~isa(R_ik, 'contSet')
                    R_ik = zonotope(R_ik);
                end
                % plot
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
            end
        end

        % plot the identification data
        y1 = squeeze(y(dim(1),k,:));
        y2 = squeeze(y(dim(2),k,:));

        if isfield(plot_settings,'s_val')
            % plot validation trajectories in gray
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

        % axis labels
        xlabel(sprintf('y_%d', dim(1)));
        ylabel(sprintf('y_%d', dim(2)));

        % set title or save as pdf
        if ~isfield(plot_settings,'path_fig')
            title("Reachable Output Sets at k=" + k);
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
