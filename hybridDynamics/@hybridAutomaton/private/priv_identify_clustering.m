function [clus,sys,data,noise] = priv_identify_clustering(traj,options)
% priv_identify_clustering - This function implements the clustering 
%   approach from Section 3.1 in [1], where the feature space used for 
%   clustering is defined by the entries of the system matrix for the 
%   dynamics identified from a local dataset constructed form the previous 
%   and next points in time.              
%
% Syntax:
%   [clus,dyn,traj,noise] = priv_identify_clustering(data,options)
%
% Inputs:
%   -traj:  trajectory object storing the trajectories used for training, 
%   -options: a structure containing user-defined algorithm parameters
%
%       .thres:  threshold for the improvement of the prediction error used
%                to stop increasing the number of clusters (Equation (8) 
%                in [1]).
%       .ARX:    order for the AutoRegressive model with eXogenous 
%                inputs (ARX).
%
% Outputs:
%   -clus:     array containing the cluster assignments for all data points
%   -sys:      cell-array storing the system dynamics for each cluster
%   -traj:     cell-array storing the modified training trajectories
%   -noise:    array containing the noise estimate for each dimension
%
% References:
%   [1] N. Kochdumper and et al. "Robust Identification of Hybrid Automata 
%       from Noisy Data", HSCC 2025
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton/identify

% Authors:       Niklas Kochdumper
% Written:       05-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % bring data to the correct format
    data = aux_initialization(traj,options);

    % identify system matrices for different sizes of local datasets
    upperBound = 10;

    Aall = cell(upperBound,1); best = cell(upperBound,1); 
    index = cell(upperBound,1); noise = cell(upperBound,1);

    for i = 1:upperBound    

        % identify system matrix
        [Aall{i},best{i},index{i},noise{i}] = aux_systemMatrixStencil(data,i);

        % update noise estimate
        if i > 1
            noise{i} = max([noise{i},noise{i-1}],[],2);
        end
    end

    % loop over different number of clusters
    for cl = 1:10

        % loop over different window sizes
        for i = 1:length(Aall)

            % loop over different ranks of the SVD
            for j = 1:max(unique(best{i}))
                
                % construct system matrix for this rank
                A = zeros(length(index{i}{1}),size(Aall{i},2));
    
                for k = 1:size(A,2)
                    A(:,k) = Aall{i}(index{i}{best{i}(j,k)},k);
                end
    
                % cluster the points based on the identified system matrices
                try
                    clus = kmeans(A',cl,'Distance','cosine','Replicates',10);
                catch
                    clus = kmeans(A',cl,'Distance','sqeuclidean','Replicates',10);
                end
        
                % compute the error for clustering
                [err(i,j,cl),sys{i,j,cl}] = aux_errorClustering(clus,data,noise{i});
                clus_{i,j,cl} = clus; 
            end
        end

        % check for convergence (see Equation (8) in [1])
        err(err == 0) = inf;

        if cl > 1
            e = squeeze(min(min(err,[],1),[],2));
            delta = (e(end-1) - e(end))/(0.5*(e(end) + e(end-1)));

            if delta <= options.thres
                cl = cl - 1; break;
            end
        end
    end

    % select the best clustering
    err = squeeze(err(:,:,cl));
    [row,col] = find(err == min(min(err)),1);

    clus = clus_{row,col,cl};
    sys = sys{row,col,cl};

    noise = noise{row};
end


% Auxiliary functions -----------------------------------------------------

function data = aux_initialization(traj,options)
% bring the data to the correct format
    
    % transform trajectory to data struct
    data = cell(length(traj),1);
    for i = 1:length(traj)
        data{i}.x = traj(i).x';
        data{i}.t = traj(i).t';
        if ~isempty(traj(i).u)            
            data{i}.u = traj(i).u';
        end
    end
    
    % remove duplicate points from trajectories
    data = aux_removeDuplicatePoints(data);

    % convert the data to discrete time
    data = aux_uniformSampledData(data);

    % add derivatives as additional states for ARX models
    if options.ARX > 1
        data = aux_dataFormatARX(data,options);
    end

    % compute the derivatives for all trajectories
    for i = 1:length(data)
        data{i}.dx = aux_derivative(data{i}.t,data{i}.x);
    end
end

function [A,best,index,noise] = aux_systemMatrixStencil(trajData,stencil)
% identify a system matrix A for each point by moving a window of fixed
% size along the trajectories and identify A using all points within the
% window

    % initialization
    n = size(trajData{1}.x,2);
    m = 0;

    if isfield(trajData{1},'u')
        m = size(trajData{1}.u,2);
    end

    a = (n+m+1)*n;
    index = cellfun(@(x) (x-1)*a+(1:a),num2cell(1:a),'UniformOutput',false);

    stencil = stencil*(n+m+1);
    wInit = ceil(stencil/2);

    A = []; noise = []; best = [];

    % loop over all trajectories
    for i = 1:length(trajData)

        Atraj = []; err = []; noi = []; bestTmp = [];
        w = min(wInit,floor((size(trajData{i}.x,1)-1)/2));

        % loop over all points of the trajectory (top of Figure 2 in [1])
        for j = w+1:size(trajData{i}.x,1) - w

            % collect data points within the window
            data = trajData{i}; ind = j-w:j+w;
            data.x = data.x(ind,:);
            data.dx = data.dx(ind,:);
            data.t = data.t(ind);

            if isfield(data,'u')
                data.u = data.u(ind,:);
            else
                data.u = [];
            end

            data.dx = aux_derivative(data.t,data.x);

            % identify system matrix
            [A_,err_,noi_,best_] = aux_identifySystemDynamics(data);

            % store system matrix and error
            Atraj = [Atraj,reshape(A_,[n*(n+m+1)^2,1])];
            err = [err,err_];
            noi = [noi,noi_'];
            bestTmp = [bestTmp,best_];
        end

        % append system matrices by values for points at start and end
        Atraj = [repmat(Atraj(:,1),[1,w]),Atraj,repmat(Atraj(:,end),[1,w])];
        err = [repmat(err(1),[1,w]),err,repmat(err(end),[1,w])];
        noi = [repmat(noi(:,1),[1,w]),noi,repmat(noi(:,end),[1,w])];
        bestTmp = [repmat(bestTmp(:,1),[1,w]),bestTmp,repmat(bestTmp(:,end),[1,w])];

        % for each data point, select system matrix resulting in the lowest 
        % error within the given window (bottom of Figure 2 in [1])
        Atraj_ = Atraj; noi_ = noi; best_ = bestTmp;

        for j = w+1:size(trajData{i}.x,1) - w
            window = j-w:j+w;
            [~,ind] = min(err(window));
            Atraj(:,j) = Atraj_(:,window(ind));
            noi(:,j) = noi_(:,window(ind));
            bestTmp(:,j) = best_(:,window(ind));
        end

        % append with matrices for other trajectories
        A = [A,Atraj]; best = [best,bestTmp];
        noise = [noise,noi(:,2*w+1:size(noi,2)-2*w-1)];
    end

    % remove outliers from the noise estimate
    ind = 1:size(noise,2);

    for i = 1:size(noise,1)
        [tmp,ind1] = sort(noise(i,:),'descend');
        m = floor(length(tmp)/2);
        ind2 = find(tmp(1:m-1) > 2*tmp(2:m));
        ind = setdiff(ind,ind1(ind2));
    end

    noise = max(noise(:,ind),[],2);
end

function [A,errBest,noise,best] = aux_identifySystemDynamics(data)
% robustly identify a linear system that fits the data best
    
    % apply dynamic mode decomposition (DMD) for all ranks of the SVD
    Alist = aux_dynamicModeDecomposition(data.x',data.dx',data.u');

    % select the matrix that best fits the data
    errBest = inf; best = zeros(length(Alist),1);

    for i = 1:length(Alist)

        % compute error
        if ~isempty(Alist{i})
            [err,noise_] = aux_computeError(Alist{i},data);
        else
            Alist{i} = zeros(size(Alist{1}));
            err = inf;
        end

        % select matrix with the lowest error
        if err < errBest
            errBest = err; best(i) = i; noise = noise_;
        else
            best(i) = best(i-1);
        end
    end

    A = horzcat(Alist{:});
end

function A = aux_dynamicModeDecomposition(X1,X2,U)
% compute the matrix X2 = A*X1 that best fits the data

    % extend matrices by entries for constant inputs and inputs
    X1 = [X1; ones(1,size(X1,2))];

    if ~isempty(U)
        X1 = [X1; U];
    end

    % singular value decomposition
    [V,S,W] = svd(X1,'econ');
    
    % construct matrices with different rank (see Equation (6) in [1])
    rankMax = sum(diag(S) > 0);
    A = cell(length(diag(S)),1);

    for rank = 1:rankMax

        % reduce rank by removing the smallest singular values
        if ~isempty(rank) && rank < size(S,1)
            V_ = V(:,1:rank); S_ = S(1:rank,1:rank); W_ = W(:,1:rank);
        else
            V_ = V; S_ = S; W_ = W;
        end
    
        % compute resulting system matrix A
        A{rank} = X2*W_*diag(1./diag(S_))*V_';
    end
end

function [err,noise] = aux_computeError(Aall,data)
% compute the error between the system approximation and the real data
% according to Equation (7) in [1]

    % construct linear system object
    n = size(Aall,1);
    A = Aall(:,1:n); c = Aall(:,n+1); B = Aall(:,n+2:end); 
    dt = data.t(2) - data.t(1);

    % loop over all time steps and integrate using Runge-Kutta-4
    x = zeros(size(data.x))';
    x(:,1) = data.x(1,:)';

    if isempty(data.u)
        for i = 1:size(x,2)-1
            k1 = A*x(:,i) + c;
            k2 = A*(x(:,i) + dt*k1/2) + c;
            k3 = A*(x(:,i) + dt*k2/2) + c;
            k4 = A*(x(:,i) + dt*k3) + c;
            x(:,i+1) = x(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        end
    else
        for i = 1:size(x,2)-1
            temp = B*data.u(i,:)' + c;

            k1 = A*x(:,i) + temp;
            k2 = A*(x(:,i) + dt*k1/2) + temp;
            k3 = A*(x(:,i) + dt*k2/2) + temp;
            k4 = A*(x(:,i) + dt*k3) + temp;
            x(:,i+1) = x(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        end
    end

    % compute error
    err = sum(sqrt(sum((x'-data.x).^2,2)));

    % compute estimate for the noise
    delta = max(abs(x(:,2:end)-x(:,1:end-1)),[],2);
    errMax = max(abs(x'-data.x),[],1);

    noise = max([errMax;delta'],[],1);
end

function [err,dyn] = aux_errorClustering(clus,trajData,noise)
% compute the averates error of the current data clustering

    % check if the system has inputs
    inputs = false;

    if isfield(trajData{1},'u')
        inputs = true;
    end

    % identify system dynamics for all clusters
    A = cell(max(clus),1); B = cell(length(A),1); c = cell(length(A),1);
    dyn = cell(length(A),1); Aall = cell(length(A),1);

    for cl = 1:length(A)
        sys = aux_constructSystemDynamics(clus,cl,trajData,noise);
        A{cl} = sys.A; c{cl} = sys.c; Aall{cl} = [sys.A,sys.c];
        if inputs
            B{cl} = sys.B; Aall{cl} = [Aall{cl},sys.B];
        end
        dyn{cl} = sys;
    end
    
    % compute error
    err = aux_trajectoryError(trajData,clus,A,B,c,inputs,noise);
end

function err = aux_trajectoryError(trajData,clus,A,B,c,inputs,noise)
% compute the error for the current cluster assignment
    
    err = 0; cnt = 1;

    % compute the error for all trajectories
    dt = trajData{1}.t(2) - trajData{1}.t(1);

    for j = 1:length(trajData)

        clus_ = clus(cnt:cnt+size(trajData{j}.x,1)-1);

        % loop over all time steps and integrate using Runge-Kutta-4
        x = zeros(size(trajData{j}.x))';
        x(:,1) = trajData{j}.x(1,:)';
    
        if ~inputs
            for i = 1:size(x,2)-1

                % check if a jump is present
                if all(abs(trajData{j}.x(i,:)' - trajData{j}.x(i+1,:)') < 2*noise)

                    A_ = A{clus_(i)}; c_ = c{clus_(i)};
    
                    k1 = A_*x(:,i) + c_;
                    k2 = A_*(x(:,i) + dt*k1/2) + c_;
                    k3 = A_*(x(:,i) + dt*k2/2) + c_;
                    k4 = A_*(x(:,i) + dt*k3) + c_;
    
                    x(:,i+1) = x(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);

                else
                    x(:,i+1) = trajData{j}.x(i+1,:)';
                end
            end
        else
            for i = 1:size(x,2)-1

                % check if a jump is present
                if all(abs(trajData{j}.x(i,:)' - trajData{j}.x(i+1,:)') < 2*noise)

                    A_ = A{clus_(i)}; c_ = c{clus_(i)}; B_ = B{clus_(i)};
                    temp = B_*trajData{j}.u(i,:)' + c_;
    
                    k1 = A_*x(:,i) + temp;
                    k2 = A_*(x(:,i) + dt*k1/2) + temp;
                    k3 = A_*(x(:,i) + dt*k2/2) + temp;
                    k4 = A_*(x(:,i) + dt*k3) + temp;
    
                    x(:,i+1) = x(:,i) + dt/6*(k1 + 2*k2 + 2*k3 + k4);

                else
                    x(:,i+1) = trajData{j}.x(i+1,:)';
                end
            end
        end

        % compute error
        err = err + sum(sqrt(sum((x'-trajData{j}.x).^2,2)));

        cnt = cnt + size(trajData{j}.x,1);
    end

    err = err/(cnt-1);
end

function [sys,data] = aux_constructSystemDynamics(clus,index,trajData,noise)
% construct the system dynamcis for the given cluster

    % get all trajectory segments belonging to the current cluster
    data = {}; cnt = 1;

    for i = 1:length(trajData)

        ind1 = find(clus(cnt:cnt+size(trajData{i}.x,1)-1) == index);

        if ~isempty(ind1)

            % find points where the cluster changes along the trajecotry
            ind2 = find(diff(ind1) > 1);
            ind2 = [0;ind2;length(ind1)];
            
            % loop over all changepoints
            for j = 1:length(ind2)-1

                ind = ind1(ind2(j)+1:ind2(j+1));
    
                if length(ind) == 1
                    if ind > 1
                        ind = [ind-1,ind];
                    end
                    if ind(end) < size(trajData{i}.x,1)
                        ind = [ind,ind(end)+1];
                    end
                end

                data{end+1} = aux_subsetStruct(trajData{i},ind);
            end
        end

        cnt = cnt + size(trajData{i}.x,1);
    end

    % split trajectories at jumps
    len = length(data);

    for i = 1:len

        % find jumps in the data
        ind = find(any(abs(data{i}.x(1:end-1,:)' - ...
                                    data{i}.x(2:end,:)') >= 2*noise,1));

        if ~isempty(ind)

            % loop over all segments between jumps
            ind = [0,ind,length(data{i}.t)];
    
            for j = 2:length(ind)-1
                if ind(j+1) ~= ind(j)+1
                    data{end+1} = aux_subsetStruct(data{i},ind(j)+1:ind(j+1));
                    data{end}.dx = aux_derivative(data{end}.t,data{end}.x);
                end
            end
    
            % consider last segment
            if ind(1)+1 ~= ind(2)
                data{i} = aux_subsetStruct(data{i},ind(1)+1:ind(2));
                data{i}.dx = aux_derivative(data{i}.t,data{i}.x);
            else
                data{i} = [];
            end
        end
    end

    data = data(~cellfun('isempty',data));

    % identify system dynamics
    if ~isempty(data)
        sys = aux_identifyLinearSys(data);
    else
        sys = [];
    end
end

function obj = aux_subsetStruct(obj,ind)
% return the subset of the struct specified by the given indices

    obj.x = obj.x(ind,:);
    obj.t = obj.t(ind);
    obj.dx = obj.dx(ind,:);

    if isfield(obj,'u')
        obj.u = obj.u(ind,:);
    end
end

function sys = aux_identifyLinearSys(data)
% identify a linear system best fitting the data

    % split the data into single data points
    points = aux_getDataPoints(data);

    % bring inputs to the correct format
    for i = 1:length(data)
        if ~isfield(data{i},'u')
            data{i}.u = [];
        end
    end

    % apply dynamic mode decomposition (DMD) for all ranks of the SVD
    Alist = aux_dynamicModeDecomposition(points.x',points.dx',points.u');

    % select the matrix that best fits the data
    errBest = inf;

    for i = 1:length(Alist)
        err = 0;
        if ~isempty(Alist{i})
            for j = 1:length(data)
                errTemp = aux_computeError(Alist{i},data{j});
                err = err + errTemp;
            end
            if err < errBest
                Aall = Alist{i};
                errBest = err;
            end
        end
    end 

    % construct linear system object
    n = size(Aall,1);
    A = Aall(:,1:n); c = Aall(:,n+1); B = Aall(:,n+2:end); 
    
    sys = linearSys(A,B,c);
end

function trajData = aux_dataFormatARX(trajData,options)
% bring data to the format required for Autoregressive Models with
% eXogenous inputs (ARX) by adding higher-order derivatives as additional
% states

    n = size(trajData{1}.x,2);

    % compute an estimate for the noise on the original data
    for i = 1:length(trajData)
        trajData{i}.dx = aux_derivative(trajData{i}.t,trajData{i}.x);
    end

    [~,~,~,noise] = aux_systemMatrixStencil(trajData,1);

    % loop over derivatives of different order that are added as states
    for i = 1:(options.ARX - 1)
        for j = 1:length(trajData)

            % compute derivative 
            dx = aux_derivative(trajData{j}.t,trajData{j}.x(:,end-n+1:end));

            % improve derivative estimate using a moving average filter
            for l = 5:30

                % terminate if derivative estimate is clean enough
                ddx = aux_derivative(trajData{j}.t,dx);
                tmp = trajData(j); tmp{1}.x = dx; tmp{1}.dx = ddx;

                [~,~,~,noise_] = aux_systemMatrixStencil(tmp,1);

                if all(noise_ <= 2*noise)
                    break;
                end

                % apply a moving average filter
                for k = 1:n
                    tmp = [dx(1,k)*ones(l,1);dx(:,k);dx(end,k)*ones(l,1)];
                    tmp = filter((1/l)*ones(1,l),1,tmp);
                    dx(:,k) = tmp(l+1:end-l);
                end
            end

            % add derivative as a new state
            trajData{j}.x = [trajData{j}.x,dx];
        end
    end
end

function points = aux_getDataPoints(trajData)
% transform the data into a list of data points

    points.x = [];
    points.dx = [];
    points.u = [];
    points.xNext = [];
    points.dt = [];

    for i = 1:length(trajData)

        m = size(trajData{i}.x,1)-1;

        points.x = [points.x; trajData{i}.x(1:end-1,:)];
        points.dx = [points.dx; trajData{i}.dx(1:end-1,:)];
        points.xNext = [points.xNext; trajData{i}.x(2:end,:)];
        
        if isfield(trajData{i},'u')
            points.u = [points.u; trajData{i}.u(1:m,:)];
        end
    end
end

function data = aux_uniformSampledData(data)
% convert the measured traces to uniform sampled traces

    % select time step size
    t = [];

    for i = 1:length(data)
        t = [t; diff(data{i}.t)];
    end

    dt = mean(t);

    if all(abs(t-dt) < 1e-10)
        return;
    end 

    % convert all measured traces to discrete time
    for i = 1:length(data)

        % detect duplicate time points to avoid interpolating between jumps
        if size(data{i}.t,2) > size(data{i}.t,1)
            data{i}.t = data{i}.t';
        end

        if length(data{i}.t) > 1 && abs(data{i}.t(end) - data{i}.t(end-1)) < eps
            data{i}.t = data{i}.t(1:end-1);
            data{i}.x = data{i}.x(1:end-1,:);
            if isfield(data{i},'u')
                data{i}.u = data{i}.u(1:end-1,:);
            end
        end

        t = (1:floor(data{i}.t(end)/dt))*dt;

        while true
            
            % recursively remove duplicate time points
            dup = find(abs(data{i}.t(2:end) - data{i}.t(1:end-1)) < eps);
            dup = [0;dup;length(data{i}.t)];

            tmp = find(dup(1:end-1) == dup(2:end)-1);

            if isempty(tmp)
                break;
            else
                indKeep = setdiff(1:length(data{i}.t),dup(tmp+1))';
                data{i}.t = data{i}.t(indKeep);
                data{i}.x = data{i}.x(indKeep,:);
                if isfield(data{i},'u')
                    data{i}.u = data{i}.u(indKeep,:);
                end
            end

        end

        % loop over all windows between duplicate time points
        x = []; u = []; tNew = []; cnt = 1;

        for j = 1:length(dup)-1
        
            % get unique indices of current window
            indOrig = dup(j)+1:dup(j+1);
            [~,tmp] = find(t <= data{i}.t(indOrig(end)));
            indNew = cnt:tmp(end);
            cnt = tmp(end) + 1;

            [~,ind] = unique(data{i}.t(indOrig));
            ind = indOrig(ind);
    
            % interpolate state vector
            x_ = interp1(data{i}.t(ind),data{i}.x(ind,:),t(indNew),'linear','extrap');
    
            if size(x_,2) ~= size(data{i}.x,2)
                x_ = x_';
            end

            x = [x;x_];
            tNew = [tNew;t(indNew)'];
    
            % interpolate input vector
            if isfield(data{i},'u')
                if size(data{i}.u,1) < size(data{i}.x,1)
                    data{i}.u = [data{i}.u;data{i}.u(end,:)];
                end
                u_ = interp1(data{i}.t(ind),data{i}.u(ind,:),t(indNew),'nearest','extrap');
                if size(u_,2) ~= size(data{i}.u,2)
                    u_ = u_';
                end
                u = [u;u_];
            end
        end

        data{i}.x = x;
        data{i}.t = tNew;

        if ~isempty(u)
            data{i}.u = u;
        end
    end
end

function trajData = aux_removeDuplicatePoints(trajData)
% remove duplicate points from the trajectories

    for i= 1:length(trajData)
        ind = 1;
        for j = 1:size(trajData{i}.x,1)-1
            if norm(trajData{i}.x(j,:)-trajData{i}.x(j+1,:)) > eps
                ind = [ind, j+1];
            end
        end
        trajData{i}.x = trajData{i}.x(ind,:);
        trajData{i}.t = trajData{i}.t(ind);
        if isfield(trajData{i},'u')
            trajData{i}.u = trajData{i}.u(ind,:);
        end
    end
end

function [dx,t] = aux_derivative(t,x)
% numerically compute the time derivative with a finite different according
% to Equation (3) in [1]
    
    % number of subintervals
    N = length(t)-1;
    
    % preallocates vector to store derivative
    dx = zeros(size(x));
    
    % approximates derivative at lower bound using forward difference
    dx(1,:) = (x(2,:)-x(1,:))/(t(2)-t(1));
    
    % approximates derivative at upper bound using backward difference
    dx(N+1,:) = (x(N+1,:)-x(N,:))/(t(N+1)-t(N));
    
    % approximates derivatives at all other nodes using central differences
    for i = 2:N
        dx(i,:) = (x(i+1,:)-x(i-1,:))/(t(i+1)-t(i-1));
    end
end

% ------------------------------ END OF CODE ------------------------------
