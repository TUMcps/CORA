function R = guardIntersect_pancake(loc,R0,guard,guardID,params,options)
% guardIntersect_pancake - implementation of the time scaling approach
%    described in [1]
%
% Syntax:
%    R = guardIntersect_pancake(loc,R0,guard,guardID,params,options)
%
% Inputs:
%    loc - location object
%    R - list of intersections between the reachable set and the guard
%    guard - guard set (class: constrained hyperplane)
%    guardID - ID of the guard set
%    params - model parameters
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - set enclosing the guard intersection
%
% References: 
%   [1] S. Bak et al. "Time-Triggered Conversion of Guards for Reachability
%       Analysis of Hybrid Automata"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       05-November-2018             
% Last update:   20-November-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % initialization
    sys = loc.contDynamics;

    % check if guard set is a constrained hyperplane
    if ~(isa(guard,'polytope') || representsa_(guard,'conHyperplane',1e-12))
        throw(CORAerror('CORA:specialError',...
            "The method 'pancake' only supports guards given as polytope objects that represent constrained hyperplanes!")); 
    end
    
    % convert hyperplane to a halfspace that represents the outside of the
    % invariant set
    c = center(R0);
    P = polytope(guard.Ae,guard.be);

    if contains_(P,c,'exact',1e-12,0,false,false)
        P = polytope(-guard.Ae,-guard.be);
    end

    % set default options for nonlinear system reachability analysis
    optionsScaled = options;
    
    if ~isa(sys,'nonlinearSys')
       optionsScaled = aux_defaultOptions(options); 
    end
    
    % create system for the time-scaled system dynamics
    [sys_,params] = aux_scaledSystem(sys,P,R0,guardID,params);

    % compute the reachable set for the time scaled system 
    R = aux_reachTimeScaled(sys_,P,R0,params,optionsScaled);
    
    % jump accross the guard set in only one time ste
    R = aux_jump(sys,P,R,params,options);

    % project the reachable set onto the hyperplane
    R = projectOnHyperplane(R,guard);

end


% Auxiliary functions -----------------------------------------------------

function [sys,params] = aux_scaledSystem(sys,P,R0,guardID,params)
% Scale the system dynamics using the distance to the hyperplane as a 
% scaling factor 

    % get maximum distance of initial set ot hyperplane
    maxDist = supremum(interval(P.A * R0 + (-P.b)));
    params.paramInt = maxDist;

    % define scaling function
    g = @(x,p) (P.A * x - P.b)./p;

    % get system dynamics
    n = sys.nrOfDims;
    m = sys.nrOfInputs;
    
    if isa(sys,'linearSys')
       f = @(x,u) aux_dynamicsLinSys(x,u,sys); 
    else
       f = sys.mFile; 
    end

    % time scaled system dynamics
    F = @(x,u,p) g(x,p) * f(x,u);

    % create symbolic variables
    xSym = sym('x',[n,1]);
    uSym = sym('u',[m,1]);
    pSym = sym('p',1);

    % create file path
    name = ['generated_',sys.name,'_',num2str(guardID),'_timeScaled'];
    foldername = [CORAROOT filesep 'models' filesep 'auxiliary'];
    path = [foldername filesep name];

    % create file for time scaled dynamics
    func = F(xSym,uSym,pSym);
    if ~isfolder(foldername)
        mkdir([CORAROOT filesep 'models'],'auxiliary');
    end
    funcHandle = matlabFunction(func,'File',path,'Vars',{xSym,uSym,pSym});

    % remove and add path so that file can be found in input argument check
    % of nonlinParamSys constructor called below
    fullpath = genpath(foldername);
    warOrig = warning;
    warning('off','all');
    rmpath(fullpath);
    warning(warOrig);
    addpath(fullpath);

    % create time scaled system
    sys = nonlinParamSys(funcHandle,n,m,1);
end

function Rfin = aux_reachTimeScaled(sys,P,R0,params,options)
% Compute the reachable set of the scaled system such that the final
% reachable set until the scaled reachable set is very close to the
% hyperplane

    % adapt options
    spec = specification(P,'unsafeSet');
    params.R0 = R0;
    if isfield(options,'maxError')
       options.maxError = inf*ones(size(options.maxError));
    end

    % prevent validateOptions error
    if isa(sys,'nonlinParamSys') && ~isfield(options,'intermediateTerms')
        options.intermediateTerms = 4;
    end
    
    % compute reachable set until 
    R = reach(sys,params,options,spec);
    
    % get final reachable set
    Rfin = R.timePoint.set{end};
end

function Rcont = aux_jump(sys,P,R0,params,options)
% compute the reachable set in such a way that the reachable set jumps in 
% only one time step accross the hyperplane    

    params.R0 = R0;
    timeStep = options.timeStep;
    timeStep_ = timeStep;
    
    % compute reachable set
    options.timeStep = timeStep;
    params.tFinal = timeStep;
    params.tStart = 0;

    R = reach(sys,params,options);

    % check if located inside the invariant
    dist_ = supportFunc_(R.timePoint.set{end},P.A','upper','interval',8,1e-3) - P.b;
    
    if dist_ < 0 
    % guard set crossed -> reduce time step size to get smaller set
        
        Rcont = R.timeInterval.set{end};
        distMin = supportFunc_(R0,P.A','lower','interval',8,1e-3) - P.b;
        lb = 0; ub = timeStep;
        
        for i = 1:10
           
            % update time step
            timeStep = (ub-lb)/2;
            options.timeStep = timeStep;
            params.tFinal = timeStep;
            
            % compute reachable set
            R = reach(sys,params,options);

            % check if located inside the invariant
            dist = supportFunc_(R.timePoint.set{end},P.A','upper','interval',8,1e-3) - P.b;
            
            if dist < 0
                Rcont = R.timeInterval.set{end};
                ub = timeStep;
                if abs(dist) <= distMin
                    break;
                end
            else
                lb = timeStep;
            end
        end
        
    else
    % guard set not crossed -> increase time interval    
        
         while true
       
            % update time step
            timeStep = timeStep + timeStep_;
            options.timeStep = timeStep;
            params.tFinal = timeStep;
            
            % compute reachable set
            R = reach(sys,params,options);

            % check if located inside the invariant
            dist = supportFunc_(R.timePoint.set{end},P.A','upper','interval',8,1e-3) - P.b;

            if dist < 0
                Rcont = R.timeInterval.set{end};
                break;
            elseif dist > dist_
                throw(CORAerror('CORA:specialError','Pancake approach failed!'));
            else
                dist_ = dist;
            end
         end 
    end
end

function options = aux_defaultOptions(options)
% set default options for nonlinear system reachability analysis (required
% if the continuous dynamics of the hybrid automaton is linear)

    % define options and default values
    opts = {'alg','tensorOrder','errorOrder','intermediateOrder', ...
            'zonotopeOrder','taylorTerms','timeStep','intermediateTerms'};
    defVal = {'lin', 3, 5, 50, 50, 10, 0.01, 4};
    
    % parse options
    list = fields(options);
    
    for i = 1:length(list)
        if ~ismember(list{i},opts)
            options = rmfield(options,list{i});
        end
    end
    
    for i = 1:length(opts)
        if ~isfield(options,opts{i})
            options.(opts{i}) = defVal{i};
        end
    end
end

function f = aux_dynamicsLinSys(x,u,sys)
% dynamic function of a linear system
    f = sys.A * x + sys.B * u + sys.c;
end
    
% ------------------------------ END OF CODE ------------------------------
