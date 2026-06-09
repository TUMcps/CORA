function traj = priv_simulateFalsify(sys,params,options)
% priv_simulateFalsify - performs several simulations of extreme system
%   behaviour using falsification algorithms
%
% Syntax:
%    res = priv_simulateFalsify(sys,params,options)
%
% Inputs:
%    sys - contDynamics object
%    params - model parameters
%    options - settings for extreme simulation via falsification
%
% Outputs:
%    res - object of class trajectory storing time and states of the 
%          simulated trajectories.

% Authors:       Niklas Kochdumper
% Written:       11-December-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    if isfield(params,'uTrans')
        params.u = params.uTrans;
    elseif isfield(params,'uTransVec')
        params.u = params.uTransVec;
    end

    list = {'R0','U','u','tFinal','tStart'};
    params = aux_removeRedundantParams(params,list);

    % split points half half between final time and total time horizon
    if mod(options.points,2) == 1
        pointsFinal = (options.points - 1)/2;
        pointsAll = options.points - pointsFinal;
    else
        pointsFinal = options.points/2;
        pointsAll = options.points/2;
    end

    % generate suitable directions
    n = sys.nrOfOutputs;
    
    if ~isfield(options,'dims')
        dims = 1:n;
    else
        dims = unique(options.dims);
    end

    dirFinal = aux_generateDirections(n,dims,pointsFinal);
    dirAll = aux_generateDirections(n,dims,pointsAll);

    % determine extreme trajectories via falsification
    traj1 = aux_falsify(sys,params,dirAll,'all');
    traj2 = aux_falsify(sys,params,dirFinal,'final');

    traj = [traj1;traj2];

end


% Auxiliary functions -----------------------------------------------------

function traj = aux_falsify(sys,params,dirs,type)
% determine the extreme trajectories in the specified directions via
% falsification

    traj = [];

    % loop over all directions
    for i = 1:size(dirs,2)

        offset = 100; res = true;

        while res

            % generate fake specification
            P = polytope(dirs(:,i)',100);
    
            if strcmp(type,'all')
                spec = specification(P,'safeSet');
            else
                spec = specification(P,'safeSet',interval(params.tFinal));
            end
    
            % determine extreme directory via falsification
            [res,fals] = falsify(sys,params,spec);

            offset = offset*10;
        end

        traj = [traj;fals.traj];
    end
end

function dirs = aux_generateDirections(n,dims,points)
% generate suitable directions for falsification

    % axis-aligned directions
    dirs = zeros(n,2*length(dims));
    dirs(dims,:) = [eye(length(dims)),-eye(length(dims))];

    % use a subset of the axis-aligned directions
    if points <= length(dims)

        ind = randperm(size(dirs,2));
        ind = ind(1:points);
        dirs = dirs(:,ind);

    % generate additional random directions
    else
        dirs_ = rand(n,points-size(dirs,2)) - 0.5;
        len = sqrt(sum(dirs_.^2,1));
        dirs_ = dirs_./len;
        dirs_ = [dirs,dirs_];
    end
end

function params = aux_removeRedundantParams(params,list)
% remove redundant parameters from struct

    fields = fieldnames(params);

    for i = 1:length(fields)
        if ~ismember(fields{i},list)
            params = rmfield(params,fields{i});
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
