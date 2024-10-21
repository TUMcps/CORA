function [res,S] = representsa_(SpS,type,tol,varargin)
% representsa_ - checks if a spectrahedral shadow can also be represented
%    by a different set representation, e.g., a special case
%
% Syntax:
%    res = representsa_(SpS,type,tol)
%    [res,Set] = representsa_(SpS,type,tol)
%
% Inputs:
%    SpS - spectraShadow object
%    type - other set representation or 'origin', 'point', 'hyperplane'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    Set - converted set
% 
% Example:
%    SpS = spectraShadow([-1 0]);
%    representsa_(SpS,'emptySet',0) % true
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa

% Authors:       Adrian Kulmburg
% Written:       02-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check empty object case
    if nargout == 1
        [empty,res] = representsa_emptyObject(SpS,type);
    else
        [empty,res,S] = representsa_emptyObject(SpS,type);
    end
    if empty; return; end

    switch type
        case 'emptySet'
            if ~isempty(SpS.emptySet.val)
                res = SpS.emptySet.val;
                return
            end
            
            % Compute the support function in any direction; if it returns
            % -Inf, the set is empty
            res = ~aux_feasible(SpS);
            if res
                SpS.emptySet.val = true;
                SpS.fullDim.val = false;
                SpS.bounded.val = true;
            end
            if nargout == 2 && res
                S = emptySet(dim(S));
            end

        case 'origin'
            n = dim(SpS);
            res = true;
            for i=1:n
                ei = zeros([n 1]);
                ei(i) = 1;
                if ~withinTol(supportFunc_(SpS,ei,'upper'),0,tol) || ~withinTol(supportFunc_(SpS,-ei,'upper'),0,tol)
                    res = false;
                    break
                end
            end
            if res
                SpS.emptySet.val = false;
                SpS.fullDim.val = false;
                SpS.bounded.val = true;
            end
            if nargout == 2 && res
                S = zeros(n,1);
            end

        case 'point'
            p = priv_findFeasiblePointSpectrahedron(SpS);
            p = SpS.c + SpS.G * p;

            SpS_shifted = SpS - p;
            res = SpS_shifted.representsa_('origin',tol);

            if res
                S = p;
            else
                S = [];
            end

        otherwise
            throw(CORAerror('CORA:specialError',...
                strcat('The function represents is not currently implemented for spectrahedra and type=',type)));
    end
end


% Auxiliary functions -----------------------------------------------------

function res = aux_feasible(SpS)

    [A0,Ai] = priv_getCoeffMatrices(SpS);
    m = size(Ai,2);
    
    % For some stupid reason, we also need to manually eliminate the case where
    % the coeff matrices Ai are all zeros
    % (otherwise, when telling Yalmip to compute Ai{i} * beta(i), it will
    % replace it by the zero matrix, not by whatever class objective they use)
    all_Ai_zero = true;
    for i = 1:m
        if ~all(all(~Ai{i}))
            all_Ai_zero = false;
            break
        end
    end
    if all_Ai_zero
        % The output now depends on whether or not A0 is PSD...
        if all(eig(A0)>=0)
            % If yes, then the problem is unbounded, but feasible
            res = true;
        else
            % If not, the problem is infeasible
            res = false;
        end
        return
    end
    
    
    beta = sdpvar(m,1,'full');
    
    A = A0;
    for i=1:m
        A = A + Ai{i} * beta(i);
    end
    
    constraints = A>=0;
    cost = [];
    persistent options
    if isempty(options)
        options = sdpsettings('solver','sedumi','verbose',0,'allownonconvex',0);
    end
    yalmipOptimizer = optimizer(constraints,cost,options,[],beta);
    
    try
        [~, exitflag] = yalmipOptimizer();
    catch ME
        if strcmp(ME.identifier,'MATLAB:nonExistentField')
            % Weird bug, specifically with SEDUMI, that means that no solutions
            % has been found
            exitflag = 1;
        else
            rethrow(ME);
        end
    end
    
    if exitflag == 1
        % unfeasible case
        res = false;
    else
        res = true;
    end

end

% ------------------------------ END OF CODE ------------------------------
