function R = reach(nlnARX,params,options,varargin)
% reach - computes the reachable sets of the nonlinear ARX system
%
% Syntax:
%    R = reach(nlnARX,params,options)
%    [R,res] = reach(nlnARX,params,options,spec)
%
% Inputs:
%    nlnARX - nonlinearARX system
%    params - parameter defining the reachability problem
%    options - options for the computation of the reachable set
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res - true/false whether specifications are satisfied
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearARX

% Authors:       Laura Luetzow
% Written:       08-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

spec = setDefaultValues({[]},varargin);

% system dimensions
p = nlnARX.n_p;
dim_y = nlnARX.nrOfOutputs;

% create initial state set from the initial measurements
if ~isfield(params,'R0')
    if isfield(options, 'armaxAlg') && strcmp(options.armaxAlg, 'exactAddition')
        params.R0 = getR0(nlnARX, params.Y0, "poly");
    else
        params.R0 = getR0(nlnARX, params.Y0);
    end
end

% options preprocessing
[params,options] = validateOptions(nlnARX,params,options);

% initialize cell array that stores the reachable sets
tVec = params.tStart:nlnARX.dt:params.tFinal;

% input vector
if size(params.u,2) == 1
    u = repmat(params.u, 1, length(tVec));
else 
    u = params.u;
end

% initialize cell array for the reachable output sets
N_k = length(tVec);
timePoint.set = cell(N_k,1);
for i = 1:p
    timePoint.set{i} = project(params.R0,[(i-1)*dim_y+1:i*dim_y]);
end

% use linReach of nonlinearSysDT
sys_DT = nonlinearSysDT(nlnARX);

% compute symbolic derivatives
if ~contains(func2str(nlnARX.mFile), "predict") && ~contains(func2str(nlnARX.mFile), "evaluate")
    derivatives(sys_DT,options);
end

switch options.armaxAlg
    case 'exactAddition'

        [~ ,U_stacked] = getStackedU(nlnARX, u, params.U, "poly");
    
        % loop over all reachablity steps
        Y_prev = params.R0;
        for k = p+1:N_k
            options.i = k;
    
            % compute next reachable set
            params.U = U_stacked{k};
            [Y_prev,options] = linReach(sys_DT,Y_prev,params,options);
    
            % replace ID's which were created for the lin error
            if ~isempty(U_stacked{k}.id)
                ID_toChange = Y_prev.id > max(U_stacked{k}.id, [], 'all') & ...
                    Y_prev.id <= max(U_stacked{end}.id, [], 'all');
                if sum(ID_toChange) > 0
                    ID_new = Y_prev.id;
                    ID_new(ID_toChange) = nlnARX.prev_ID+1 : nlnARX.prev_ID+sum(ID_toChange);
                    setPrevID(nlnARX, nlnARX.prev_ID+sum(ID_toChange));
                    Y_prev = replaceId(Y_prev, ID_new);
                end
            end
    
            timePoint.set{k} = project(Y_prev,[(p-1)*dim_y+1:p*dim_y]);
    
            % log information
            verboseLog(options.verbose,k,tVec(k),params.tStart,params.tFinal);
    
            % check specification
            if ~isempty(spec)
                if ~check(spec,timePoint.set{k},interval(tVec(k)))
                    timePoint.set = timePoint(1:k);
                    timePoint.time = num2cell(tVec(1:k)');
                    R = reachSet(timePoint);
                    return;
                end
            end
        end

    otherwise
    % linearize and use ARX reformulation
    
        [~ ,U_stacked] = getStackedU(nlnARX, u, params.U); 
    
        % compute symbolic derivatives
        if ~contains(func2str(nlnARX.mFile), "predict") && ~contains(func2str(nlnARX.mFile), "evaluate")
            derivatives(nlnARX,options);
        end
    
        % compute linerized parameters
        p_GO = computeGO(nlnARX, params.R0, u+center(params.U), N_k);
    
        E = [zeros(dim_y,(p-1)*dim_y) eye(dim_y)];
        L = cell(N_k,1);
    
        for k = p+1:N_k
            % compute lin error set L(k) using linearSysDT
            U_save = params.U; % original input set required for Y(k)
            params.U = U_stacked{k};
            % if ~contains(func2str(sys.mFile), "predict") && ~contains(func2str(sys.mFile), "evaluate")
            try
                Y_prev = getR0(nlnARX, timePoint.set(k-p:k-1));
                [~,~,Verror] = linReach(sys_DT,Y_prev,params,options);
                L{k} = E * Verror;
                % else
                %     L{k} = zeros(dim_y, 1);
                % end
    
                % compute output set Y(k)
                params.U = U_save;
                Y_d = p_GO.y(:,k) + p_GO.C{k} * (params.R0-center(params.R0));
                for j = 1:k
                    Y_d = Y_d + p_GO.D{k,j} * (params.U-center(params.U));
                    if j >= p+1
                        Y_d = Y_d + p_GO.E{k,j} * L{j};
                    end
                end
            catch
                Y_d = zonotope([]);
            end
    
            timePoint.set{k} = Y_d;
    
            % log information
            verboseLog(options.verbose,k,tVec(k),params.tStart,params.tFinal);
    
            % check specification
            if ~isempty(spec)
                if ~check(spec,timePoint.set{k},interval(tVec(k)))
                    timePoint.set = timePoint(1:k);
                    timePoint.time = num2cell(tVec(1:k)');
                    R = reachSet(timePoint);
                    return;
                end
            end
        end
end

% create reachable set object
timePoint.time = num2cell(tVec(1:end)');
R = reachSet(timePoint);

% log information
verboseLog(options.verbose,length(tVec),tVec(end),params.tStart,params.tFinal);

% ------------------------------ END OF CODE ------------------------------
