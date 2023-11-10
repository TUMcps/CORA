function [Rerror,options] = errorSolution_adaptive(obj,options,Vdyn,Vstat)
% errorSolution_adaptive - computes the solution due to the abstraction error
%    the number of Taylor terms is chosen according to the set size decrease
%
% Syntax:
%    Rerror = errorSolution_adaptive(obj,options,Vdyn,Vstat)
%
% Inputs:
%    obj - linearized system
%    options - options struct (for nonlinear system)
%    Vdyn - set of dynamic errors
%    Vstat - set of static errors
%
% Outputs:
%    Rerror - reachable set due to the linearization error
%    options - options struct (for nonlinear system)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Mark Wetzlinger
% Written:       24-April-2020
% Last update:   27-May-2020
%                15-June-2020 (include Vstat)
%                10-July-2020 (delete Vstat from convergence process)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if static error given
if nargin < 4 || representsa_(Vstat,'emptySet',1e-10)
    isVstat = false;
else
    isVstat = true;
end
% time step
deltat = options.timeStep;

% exponential
A = obj.A;
A_abs = abs(A);
Apower{1} = A;  
Apower_abs{1} = A_abs; 
M = eye(obj.dim);
eAabst = expm(A_abs*deltat);
%initialize Asum
AVsum = deltat * Vdyn;
RerrorInt_etanoF = sum(abs(generators(AVsum)),2);
if isVstat
    Asum = deltat * eye(obj.dim);
end

eta = 1; breakcond = false;
% loop over increasing Taylor terms
while true

    %compute powers
    M = M + Apower_abs{eta} * deltat^(eta) / factorial(eta);
    
    % compute powers
    temp = deltat^(eta+1) / factorial(eta+1) * Apower{eta};
    ApowerV = temp * Vdyn;   
    % compute sum
    AVsum = AVsum + ApowerV;
    
    % compute error set
    if isVstat
        % including static error
        Asum = Asum + temp;
    end
    
    RerrorInt_etanoF = RerrorInt_etanoF + sum(abs(generators(ApowerV)),2);
    
    % at least two sets for comparsion needed
    if eta > 1
        gainnoF = 1 - formerEdgenoF / RerrorInt_etanoF(critDimnoF);

        % break if gain too small (or same truncation order as in previous
        % iteration of the same time step reached)
        if isfield(options,'tt_err') && length(options.tt_err) == options.i
            breakcond = eta == options.tt_err(options.i);
        elseif gainnoF < options.zetaTabs || formerEdgenoF == 0 % gain < thr
            breakcond = true;
        end

        if breakcond
            % Vdyn: take former set (less generators)
            % determine error due to finite Taylor series
            W = abs(eAabst - M);
            E = interval(-W,W);
            % get error due to finite Taylor series
            F = E*Vdyn*deltat;
            
            if isVstat
                % also former Asum due to E
                eAtInt = Asum + E*deltat;
                Rerror = AVsum + F + eAtInt * Vstat;
            else
                Rerror = AVsum + F;
            end
            % save taylor order for analysis
            options.tt_err(options.i,1) = eta;
            break
        end
        formerEdgenoF = RerrorInt_etanoF(critDimnoF);
    else
        [formerEdgenoF,critDimnoF] = max(RerrorInt_etanoF);
    end
    
    % compute powers
    Apower{eta+1,1} = Apower{eta}*A;
    Apower_abs{eta+1,1} = Apower_abs{eta}*A_abs;
    % increment Taylor order
    eta = eta + 1;
    
end


end

% ------------------------------ END OF CODE ------------------------------
