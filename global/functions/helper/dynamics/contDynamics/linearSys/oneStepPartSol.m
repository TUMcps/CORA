function [partSol,expmat] = oneStepPartSol(Z,expmat,dt,varargin)
% oneStepPartSol - computes an inner-approximation or
%    outer-approximation of the particular solution due to time-varying
%    inputs over one step
%    note: set Z has to be centered at the origin!
%
% Syntax:
%    [partSol,expmat] = oneStepPartSol(Z,expmat,dt,varargin)
%
% Inputs:
%    Z - zonotope object
%    expmat - exponentialMatrix object
%    dt - time step size
%    type - (optional) 'outer' or 'inner'
%
% Outputs:
%    partSol - particular solution at time dt
%    expmat - exponentialMatrix object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       13-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default: outer-approximation
type = setDefaultValues({'outer'},varargin);
if ~any(strcmp(type,{'outer','inner'}))
    throw(CORAerror('CORA:wrongValue','fourth',...
        'has to be ''outer'' or ''inner'''));
end

% state dimension
n = length(expmat.A);

% outer-approximation or inner-approximation?
switch type

    case 'outer'

        % read out center
        c = center(Z);
        c_res = c * dt;
        
        % read out generator matrix
        G = generators(Z);
        G_res = G * dt;
        G_diag = sum(abs(G_res),2);
        
        % loop until additional terms are small enough
        eta = 1;
        while true
        
            % compute powers of A
            expmat = next_Apower(expmat,eta);

            % next matrix
            % note that expmat.Apower{eta} = A^eta / eta!
            M = expmat.Apower{eta} / (eta+1) * dt^(eta+1);

            % include additional contribution to center
            c_res = c_res + M*c;
            
            % additional term
            G_add = M*G;
            G_add_diag = sum(abs(G_add),2);
        
            % concatenate generator matrix
            G_res = [G_res, G_add];
        
            % check if floating-point precision reached
            if all( abs(G_add_diag) <= eps * abs(G_diag) )
                break
            end
        
            % add term to simplified value for convergence
            G_diag = G_diag + G_add_diag;
        
            % increment Taylor order
            eta = eta + 1;
        end
        
        % instantiate particular solution (center has to be at origin)
        partSol = zonotope(c_res,G_res);

    case 'inner'
        
        % compute particular solution
        if expmat.isAinv
            partSol = expmat.Ainv * (expm(expmat.A*dt) - eye(n)) * Z;

        else
            % loop until additional terms are small enough
            eta = 1;
            M = eye(n) * dt;
            while true
            
                % compute powers of A
                expmat = next_Apower(expmat,eta);
                
                % additional term: A^eta * dt^(eta+1) / (eta+1)!
                % note that expmat.Apower{eta} = A^eta / eta!
                M_add = expmat.Apower{eta} / (eta+1) * dt^(eta+1);
            
                % add to solution
                M = M + M_add;
            
                % check convergence (up to floating-point precision)
                if all(all( abs(M_add) <= abs(M) * eps ))
                    break
                end
            
                % increment Taylor order
                eta = eta + 1;
            end
            
            % multiply resulting matrix with zonotope
            partSol = M * Z;
        end

end

end

% ------------------------------ END OF CODE ------------------------------
