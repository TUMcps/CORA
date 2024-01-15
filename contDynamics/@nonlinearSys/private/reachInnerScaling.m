function [Rin,Rout] = reachInnerScaling(sys,params,options)
% reachInnerScaling - compute an inner-approximation of the reachable set
%                     with the algorithm in [1]
%
% Syntax:
%    [Rin,Rout] = reachInnerScaling(sys,params,options)
%
% Inputs:
%    sys - nonlinearSys object
%    params - parameters defining the reachability problem
%    options - struct containing the algorithm settings
%
% Outputs:
%    Rin - object of class reachSet storing the inner-approximation of the 
%          reachable set
%    Rout - object of class reachSet storing the outer-approximation of the
%           reachable set
%
% References:
%    [1] N. Kochdumper and M. Althoff. "Computing Non-Convex Inner-
%        Approximations of Reachable Sets for Nonlinear Continuous Systems"
%        CDC 2020
%    [2] N. Kochdumper et al. "Utilizing Dependencies to Obtain Subsets of 
%        Reachable Sets", HSCC 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachInner

% Authors:       Niklas Kochdumper
% Written:       14-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % options preprocessing
    input = false;
    if isfield(params,'U')
        input = true;
    else
        params.U = interval(zeros(sys.nrOfInputs));
    end
    
    options = validateOptions(sys,mfilename,params,options);

    % compute outer-approximation of the reachable set
    [params,options_outer] = aux_getOuterReachOptions(options);
    
    if ~input
        Rout = reach(sys,params,options_outer);
        Rout_ = Rout;
    else
        Rout_ = aux_reachWithInputs(sys,params,options_outer,options);
        if nargout > 1
            Rout = reach(sys,params,options_outer);
        end
    end

    % loop over all time steps for the inner-approximation
    set = {options.R0}; time = {options.tStart};
    
    for i = options.N:options.N:length(Rout_.timePoint.set)
       
        % Step 1: Outer-Approximation -------------------------------------
        
        % get outer approximation
        R = Rout_.timePoint.set{i};
        
        % remove additional dependent generators
        R = aux_removeAddGens(R);
        
        % order reduction for the dependent part
        R = aux_reduceOrderDep(R,options.orderInner);
        
        
        % Step 2: Boundary Enclosure --------------------------------------
        
        % loop over the 2n faces of the initial set
        n = size(R.E,1);
        B = cell(2*n,1);
        
        for j = 1:n
            
            % extract enclosure for the boundary from the outer-
            % approximation of the reachable set using the approach in [2]
            B1 = getSubset(R,j,1);
            B2 = getSubset(R,j,-1);
            
            % reduce the order of the independent part of the poly zonotope
            B{2*(j-1)+1} = aux_reduceOrderInd(B1);
            B{2*(j-1)+2} = aux_reduceOrderInd(B2);
        end
        
        
        % Step 3: Inner-Approximation -------------------------------------
        
        % initial factor domain
        I = interval(-ones(n,1),ones(n,1));
        
        % initial scaling factor
        if strcmp(options.scaleFac,'auto')
            factor = 0.99;
        else
            factor = options.scaleFac;
        end
        
        % loop until a non-empty set is found
        loop = 1;
        
        while factor > 0 && loop
        
            loop = 0;
            
            % loop over all boundary sets
            for j = 1:length(B)

                % constraint for intersection of boundary and reach. set
                conFun = aux_constraintIntersection(R,B{j});

                % compute an initial solution using nonlinear programming
                k = ceil(j/2);

                if mod(j,2) == 0
                    Itemp = aux_initialSolution(conFun,I,B{j},factor,k,-1);
                else
                    Itemp = aux_initialSolution(conFun,I,B{j},factor,k,1);
                end

                % prove that the solution is valid by contracting the the
                % resulting domain to the empty set
                I_ = contract(conFun,Itemp,options.contractor, ...
                              options.iter,options.splits);

                if isempty(I_)                          % valid solution
                   I = project(Itemp,1:n);
                elseif strcmp(options.scaleFac,'auto')  % decrease scaling
                   factor = factor - 0.01; 
                   loop = 1; break;
                else                                    % set is empty
                   I = [];
                   warning('Inner-approximation is empty!'); 
                   break;
                end
            end
        end
        
        % construct the inner-approximation of the polynomial zonotope by
        % inserting the contractor factor domain
        if ~representsa_(I,'emptySet',eps)
            R = getSubset(R,1:n,I);
            R = polyZonotope(R.c,R.G,[],R.E);
        else
            R = [];
        end
        
        if input && ~representsa_(R,'emptySet',1e-12)
            set{end+1} = project(R,1:sys.dim);
        else
            set{end+1} = R;
        end
        time{end+1} = Rout_.timePoint.time{i+1};
    end
    
    % construct reachSet object for inner-approximation
    timePoint.set = set;
    timePoint.time = time;
    
    Rin = reachSet(timePoint);
end
    

% Auxiliary functions -----------------------------------------------------

function I0 = aux_initialSolution(conFun,I,B,factor,ind,type)
    
    % construct domain
    p = size(B.E,1) + size(B.GI,2);
    I0 = cartProd(I,interval(-ones(p,1),ones(p,1)));
    
    % compute an initial solution by optimization
    if type == 1
        fun = @(x) x(ind);
    else
        fun = @(x) -x(ind);
    end
    
    conFun_ = @(x) deal(0,conFun(x));
    x0 = center(I0);
    
    opt = optimset('Display','off');
    w = warning;
    warning('off');
    
    [~,val] = fmincon(fun,x0,[],[],[],[],infimum(I0),supremum(I0), ...
                      conFun_,opt);

    warning(w);
                  
    % scale domain by factor to make the contraction easier
    val = factor*val;
    
    % check if the solution is valid by contracting the domain to the empty
    % set
    if type == 1
        I0(ind) = interval(infimum(I0(ind)),val);
    else
        I0(ind) = interval(-val,supremum(I0(ind)));
    end
end

function conFun = aux_constraintIntersection(pZ1,pZ2)
% construct generator and exponent matrix for the polynomial constraints
% G * alpha^E + b = 0 resulting from the intersection of pZ1 and pZ2

    b = pZ1.c - pZ2.c;
    G = [pZ1.G, -pZ2.G, -pZ2.GI];
    
    p_ = size(pZ2.GI,2);
    E = blkdiag(pZ1.E, pZ2.E, eye(p_));
    
    conFun = @(x) aux_constrainedFunction(x,b,G,E);
end

function val = aux_constrainedFunction(x,b,G,E)
% constraint function for the constraint G * x^E + b = 0

    val = b;

    for i = 1:size(G,2)
        temp = 1;
        for j = 1:length(x)
           temp = temp.* x(j)^E(j,i); 
        end
        val = val + G(:,i).*temp;
    end    
end

function pZ = aux_reduceOrderDep(pZ,order)
% reduce the order of the dependent part of the polynomial zonotope

    n = length(pZ.c);
    m = size(pZ.G,2);
    
    if m > n * order
        temp = polyZonotope(pZ.c,pZ.G,[],pZ.E);
        temp = reduce(temp,'girard',order+1);
        pZ = polyZonotope(temp.c,temp.G,[temp.GI,pZ.GI],temp.E);
    end
end

function pZ = aux_reduceOrderInd(pZ)
% reduce the order of the independent part of the polynomial zonotope to 1

    % get orientation by applying Principal Component Analysis to the
    % generator matrix
    [B,~,~] = svd([-pZ.G,pZ.G]);
    
    % compute interval enclosure in the transformed space
    zono = zonotope([zeros(length(pZ.c),1),pZ.GI]);
    zono = B * zonotope(interval(B'*zono));

    % construct the resulting polynomial zonotope
    pZ = polyZonotope(pZ.c,pZ.G,generators(zono),pZ.E);
end

function pZ = aux_removeAddGens(pZ)
% remove additional dependent generators for factors that do not belong to
% the initial set but were introduced by restructuring

    % check if there are additional factors that need to be removed
    n = length(pZ.c);
    
    if size(pZ.E,1) > n
        
       % determine indices of the generators that belong to add. factors
       ind = find(sum(pZ.E(n+1:end,:),1) >= 1);
       ind_ = setdiff(1:size(pZ.G,2),ind);
       
       % enclose generators belonging to add. genertors with a zonotope
       temp = polyZonotope(zeros(n,1), pZ.G(:,ind), [], pZ.E(:,ind));
       zono = zonotope(temp);
       
       % construct the resulting polynomial zonotope
       pZ = polyZonotope(pZ.c+zono.c, pZ.G(:,ind_), ...
           [pZ.GI,zono.G], pZ.E(1:n,ind_)); 
    end
end

function [params,options] = aux_getOuterReachOptions(options)
% extract params and options for outer-reachability analysis

    % copy relevant fields to the params struct
    params.R0 = options.R0;
    params.tFinal = options.tFinal;
    if isfield(options,'U')
       params.U = zonotope(options.U);
    end
    if isfield(options,'u')
       params.u = options.u; 
    end
    if isfield(options,'tStart')
       params.tStart = options.tStart;
    end
    
    % remove irrelevant fields from the options struct
    list = {'R0','U','u','tStart','tFinal','N','orderInner','algInner', ...
            'contractor','iter','splits','scaleFac','timeStepInner', ...
            'inpChanges','linAlg'};
        
    for i = 1:length(list)
       if isfield(options,list{i})
          options = rmfield(options,list{i}); 
       end
    end
    
    options.polyZono = rmfield(options.polyZono,'volApproxMethod');
end

function Rout = aux_reachWithInputs(sys,params,options_outer,options)
% compute outer-approximation of the reachable set with piecewise constant
% inputs to obtain a more accurate inner-approximation

    inputChanges = options.inpChanges + 1;
    m = sys.nrOfInputs;
    n = sys.dim;
    
    % define extended initial set
    R0 = interval(params.R0);
    for i = 1:inputChanges
       R0 = cartProd(R0,interval(params.U));
    end
    params.R0 = polyZonotope(R0);
    params = rmfield(params,'U');
    
    t = linspace(params.tStart,params.tFinal,inputChanges+1);
    Rout = []; cnt = 1;
    
    for i = 1:inputChanges
        
        % define extended system dynamics
        fun = @(x,u) [sys.mFile(x,x(n+(i-1)*m+1:n+i*m )); ...
                      zeros(inputChanges*m,1)];
        sysExt = nonlinearSys([sys.name,'Extended', ...
                              num2str(inputChanges),'_',num2str(i)],fun);
        
        % compute reachable set
        if ~isfield(params,'u') || all(all(params.u == 0))
            
            params.tStart = t(i);
            params.tFinal = t(i+1);

            Rtemp = reach(sysExt,params,options_outer);

            Rout = add(Rout,Rtemp);
            params.R0 = Rtemp.timePoint.set{end};      
            
        else
            
            t_ = t(i):options.timeStep:t(i+1);
            
            for j = 1:length(t_)-1
                
                params.tStart = t_(j);
                params.tFinal = t_(j+1);
                params.R0 = params.R0 + [zeros(n,1);repmat(params.u(:,cnt), ...
                                        [inputChanges,1])];

                Rtemp = reach(sysExt,params,options_outer);

                Rout = add(Rout,Rtemp);
                params.R0 = Rtemp.timePoint.set{end};  
                cnt = cnt + 1;
            end
        end
    end
    
    % bring reachable set object to correct format
    timePoint.set = Rout(1).timePoint.set;
    timePoint.time = Rout(1).timePoint.time;
    timeInterval.set = Rout(1).timeInterval.set;
    timeInterval.time = Rout(1).timeInterval.time;
    
    for i = 2:length(Rout)
        timePoint.set = [timePoint.set; Rout(i).timePoint.set];
        timePoint.time = [timePoint.time; Rout(i).timePoint.time];
        timeInterval.set = [timeInterval.set; Rout(i).timeInterval.set];
        timeInterval.time = [timeInterval.time; Rout(i).timeInterval.time];
    end
    
    Rout = reachSet(timePoint,timeInterval);
end

% ------------------------------ END OF CODE ------------------------------
