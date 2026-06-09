function val = robustness(spec,varargin)
% robustness - computes the robustness score of a point with respect to 
%    the specifications, where a positive robustness score means that all
%    specifications are satisfied
%
% Syntax:
%    val = robustness(spec,traj)
%    val = robustness(spec,x)
%    val = robustness(spec,x,time)
%
% Inputs:
%    spec - specification object
%    traj - simulated trajectories (class trajectory)
%    x - states of the trace (dimensions: [m,n])
%    t - times of the trace (dimensions: [m,1])
%
% Outputs:
%    val - robustness value for the point p
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Authors:       Niklas Kochdumper
% Written:       27-November-2021             
% Last update:   17-November-2025 (NK, changed interface + restructured)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % initialize robustness value
    val = Inf;

    % catch case with multiple specifications
    if size(spec,1) > 1

        for i = 1:size(spec,1)
            val_ = robustness(spec(i,1),varargin{:});
            val = min(val,val_);
        end
        
        return;
    end

    % parse input arguments
    t = [];

    if nargin == 2
        if isa(varargin{1},'trajectory')

            traj = varargin{1};
           
            % case with multipe trajectories
            if size(traj,1) > 1
                for i = 1:size(traj,1)
                    val_ = robustness(spec,traj(i,1));
                    val = min(val,val_);
                end
                
                return;

            % case with one single trajectory
            else
                if ~isempty(traj.y)
                    x = traj.y;
                else
                    x = traj.x;
                end
                t = traj.t;
            end
        else
            x = varargin{1};
        end
    else
       x = varargin{1}; t = varargin{2};
    end

    % catch case with temporal logic specification
    if strcmp(spec.type,'logic')
        if isempty(t)
            throw(CORAerror('CORA:specialError',...
                  'Time is required for temporal logic specifications."'));
        end

        val = robustness(spec.set,x',t);
        return;
    end

    % check if time is valid
    if ~representsa_(spec.time,'emptySet',eps)
        if isempty(t)
            throw(CORAerror('CORA:specialError',...
                        'Time is required for timed specifications."'));

        % interpolate to obtain trajectory at required time
        elseif ~any(contains_(spec.time,t,'exact',eps,0,false,false))
            [~,ind] = unique(t);
            x = interp1(t(ind),x(:,ind)',center(spec.time))';
            t = center(spec.time);
        end
    end
    
    % loop over all points
    for i = 1:size(x,2)

        if representsa_(spec.time,'emptySet',eps) || ...
                     contains_(spec.time,t(i),'exact',eps,0,false,false)

            % different types of specifications
            switch spec.type

                case 'invariant'
                    val_ = aux_robustnessSafeSet(spec.set,x(:,i));

                case 'unsafeSet'
                    val_ = aux_robustnessUnsafeSet(spec.set,x(:,i));

                case 'safeSet'
                    val_ = aux_robustnessSafeSet(spec.set,x(:,i));

                case 'custom'
                    throw(CORAerror('CORA:notSupported',...
                        ['Robustness computation for custom ' ...
                         'specifications is not supported!']));
            end

            % overall robustness is minimum of single specifications
            if ~isempty(t)
                val = min(val,val_);
            else
                val = [val,val_];
            end
        end
    end

    % assign output arguments
    if isempty(t)
        val = val(2:end);
    end
end


% Auxiliary functions -----------------------------------------------------

function val = aux_robustnessUnsafeSet(S,p)
% compute the robustness value of point p for an unsafe set S

    % convert set S to a polytope with normalized halfspace directions
    S_ = polytope(S);
    
    if ~S_.isHRep.val
        constraints(S_);
    end

    S_ = normalizeConstraints(S_,'A');
    C = S_.A; d = S_.b;
    
    % check if the point is inside or outside the unsafe set
    if contains(S,p)
        val = -min(abs(C*p-d));             % distance to polytope boundary
    else
        val = aux_distancePolyPoint(C,d,p);
    end
end

function val = aux_robustnessSafeSet(S,p)
% compute the robustness value of point p for a safe set S

    % convert set S to a polytope with normalized halfspace directions
    S_ = polytope(S);
    
    if ~S_.isHRep.val
        constraints(S_);
    end

    S_ = normalizeConstraints(S_,'A');
    C = S_.A; d = S_.b;
    
    % check if the point is inside or outside the safe set
    if contains(S,p)
        val = min(abs(C*p-d));              % distance to polytope boundary
    else
        val = -aux_distancePolyPoint(C,d,p);
    end
end

function d = aux_distancePolyPoint(C,d,p)
% compute the norm 1 distance between a point p and a polytope P: C*x <= d    

    % get polytope properties
    n = length(p); m = size(C,1);
    
    % check how many halfspace constraints are violated
    halfspaceVals = C*p - d;
    ind = find(halfspaceVals > 0);
    
    if length(ind) == 1
        
        % only one halfspace constraint violated -> distance to polytope is
        % equal to the distance to the halfspace constraint
        d = halfspaceVals(ind(1));
        
    elseif length(ind) == n
        
        % compute the vertex that is closest to the point by combining 
        % the violated halfspace constraints
        v = linsolve(C(ind,:),d(ind));
        d = sqrt(sum((v-p).^2));
        
    else
        % set-up linear program to minimize the norm 1 distance: 
        % min ||p - x||_1 s.t. C*x <= d
        problem.f = [zeros(n,1); ones(2*n,1)];
        problem.Aeq = [eye(n) eye(n) -eye(n)]; problem.beq = p;
        problem.Aineq = [C zeros(m,2*n); zeros(2*n,n) -eye(2*n)];
        problem.bineq = [d; zeros(2*n,1)];
        problem.lb = [];
        problem.ub = [];

        % solve linear program
        [~,d] = CORAlinprog(problem);
    end
end

% ------------------------------ END OF CODE ------------------------------
