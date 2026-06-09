function [res,fals] = falsifyMonteCarlo(fun,params,options,spec)
% falsifyMonteCarlo - find a falsifying trajectory for a black-box system
%    using Monte-Carlo sampling
%
% Syntax:
%    [res,fals] = falsifyMonteCarlo(fun,params,options,spec)
%
% Inputs:
%    fun - function handle [t,x] = f(x0,u) for the black box simulation
%          function, where x0 is the initial state, u are the system 
%          inputs, and x and t are the states and time points of the
%          resulting trajectory
%    params - parameter defining the reachability problem
%       .R0: initial set (class contSet)
%       .U:  input set (class contSet)
%    options - options for falsification
%       .nrConstInp: number of piecewise-constant input segments 
%                    (default: [], number is determined automatically)
%       .maxTime:    maximum computation time allocated for falsification 
%                    in seconds (default: 600)
%    spec - object of class specification (reach-avoid) or stl
%
% Outputs:
%    res - true/false whether falsification was successfull
%    fals - struct containing falsifying trajectory
%           .x0   ... point from initial set
%           .u    ... piecewise-constant input values
%           .tu   ... switching times of .u
%           .traj ... object of class trajectory
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/falsify

% Authors:       Niklas Kochdumper
% Written:       24-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % loop until a falsifying trajectory is found
    tComp = 0; r = Inf; res = false; fals = []; cnt = 1;

    while true
        
        clock = tic();

        % sample random initial state and random input
        x0 = randPoint(params.R0);
        u = randPoint(params.U,options.nrConstInp);

        % combine with time-varying input
        if isfield(params,'u')

            t = linspace(0,1,size(params.u,2)+1);

            if size(params.u,2) > size(u,2)
                tu = linspace(t(1),t(end),size(params.u,2)+1);
                t_ = linspace(t(1),t(end),size(u,2)+1);
                u_ = interp1(t_,[u,u(:,end)]',tu(1:end-1),'previous')';
                u = parmas.u + u_;
            else
                tu = linspace(t(1),t(end),size(u,2)+1);
                t_ = linspace(t(1),t(end),size(params.u,2)+1);
                u_ = interp1(t_,[params.u,params.u(:,end)]', ...
                                              tu(1:end-1),'previous')';
                u = u + u_;
            end
        end

        % simulate the system
        [t,x] = fun(x0,u);

        % compute the robustness
        traj = trajectory([],x,[],t);
        r_ = robustness(spec,traj);

        % update best robustness
        if r_ < r
            r = r_; fals.x0 = x0; fals.u = u; fals.traj = traj;
            tu = linspace(t(1),t(end),size(u,2)+1); fals.tu = tu(1:end-1); 

            if r < 0
               res = true; break;
            end

            if options.verbose
               disp(['Lowest robustness: ',num2str(r), ...
                                        '  (',num2str(cnt),' samples)']); 
            end
        end

        % terminate if maximum time is exceeded
        tComp = tComp + toc(clock);

        if tComp > options.maxTime
            break;
        end

        cnt = cnt + 1;
    end

    % display final results
    if options.verbose
        if res
            disp('Falsification successfull');
        elseif tComp > options.maxTime
            disp('Stopping because time exceeded options.maxTime');
        end
    end

% ------------------------------ END OF CODE ------------------------------
