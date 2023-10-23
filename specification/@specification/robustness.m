function val = robustness(spec,p,varargin)
% robustness - computes the robustness score of a point with respect to 
%    the specifications, where a positive robustness score means that all
%    specifications are satisfied
%
% Syntax:
%    val = robustness(spec,p)
%    val = robustness(spec,p,time)
%
% Inputs:
%    spec - specification object
%    p - point represented by a vector
%    time - scalar representing the current time
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % initialize robustness value
    val = Inf;

    % parse input arguments
    time = [];
    
    if nargin > 2 && ~isempty(varargin{1})
       time = varargin{1}; 
    end
    
    % check if multiple points are provided
    if size(p,2) > 1
        
        % compute robustness for all points
        val = zeros(1,size(p,2));
        
        for i = 1:size(p,2)
            val(i) = robustness(spec,p(:,i),varargin{:}); 
        end
        
    else

        % loop over all specifications
        for i = 1:size(spec,1)

            % check if time frames overlap
            if isempty(time) && ~representsa_(spec(i,1).time,'emptySet',eps)
                throw(CORAerror('CORA:specialError',...
                    'Timed specifications require a time interval."'));
            end

            if representsa_(spec(i,1).time,'emptySet',eps) ...
                || contains(spec(i,1).time,time)

                % different types of specifications
                switch spec(i,1).type

                    case 'invariant'
                        val_ = aux_robustnessSafeSet(spec(i,1).set,p);

                    case 'unsafeSet'
                        val_ = aux_robustnessUnsafeSet(spec(i,1).set,p);

                    case 'safeSet'
                        val_ = aux_robustnessSafeSet(spec(i,1).set,p);

                    case 'custom'
                        throw(CORAerror('CORA:notSupported',...
                            ['Robustness computation for custom ' ...
                             'specifications is not supported!']));

                    case 'logic'
                        throw(CORAerror('CORA:notSupported',...
                            ['Robustness computation for logic ' ...
                             'specifications is not supported!']));
                end

                % overall robustness is minimum of single specifications
                val = min(val,val_);
            end
        end
    end
end


% Auxiliary functions -----------------------------------------------------

function val = aux_robustnessUnsafeSet(S,p)
% compute the robustness value of point p for an unsafe set S

    % convert set S to a polytope with normalized halfspace directions
    S_ = normalizeConstraints(polytope(S),'A');
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
    S_ = normalizeConstraints(polytope(S),'A');
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
    temp = C*p - d;
    ind = find(temp > 0);
    
    if length(ind) == 1
        
        % only one halfspace constraint violated -> distance to polytope is
        % equal to the distance to the halfspace constraint
        d = temp(ind(1));
        
    elseif length(ind) == n
        
        % compute the vertex that is closest to the point by combining 
        % the violated halfspace constraints
        v = linsolve(C(ind,:),d(ind));
        d = sqrt(sum((v-p).^2));
        
    else
        % set-up linear program to minimize the norm 1 distance: 
        % min ||p - x||_1 s.t. C*x <= d
        f = [zeros(n,1); ones(2*n,1)];
        Aeq = [eye(n) eye(n) -eye(n)]; beq = p;
        A = [C zeros(m,2*n); zeros(2*n,n) -eye(2*n)];
        b = [d; zeros(2*n,1)];

        % solve linear program
        options = optimoptions('linprog','Display','off');
        [~,d] = linprog(f,A,b,Aeq,beq,[],[],options);
    end
end

% ------------------------------ END OF CODE ------------------------------
