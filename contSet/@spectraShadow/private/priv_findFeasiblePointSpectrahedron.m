function p = priv_findFeasiblePointSpectrahedron(SpS)
% priv_findFeasiblePointSpectrahedron - returns a feasible point of the 
%    'base' spectrahedron of SpS. Important: This is NOT a point in SpS,
%    but rather a point in the 'base' spectrahedron of SpS!
%
% Syntax:
%    p = priv_findFeasiblePointSpectrahedron(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    p - point in the 'base' spectrahedron of SpS
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Adrian Kulmburg
% Written:       26-August-2023 
% Last update:   ---    
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    [A0,Ai] = priv_getCoeffMatrices(SpS);
    m = size(Ai,2);

    % For some stupid reason, we also need to manually eliminate the case,
    % where the coeff matrices Ai are all zeros
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
        % We know that the problem has to be feasible, so here we have to
        % assume that any point will do the job; we take the origin
        p = zeros([dim(SpS) 1]);
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
    [p, exitflag] = yalmipOptimizer();

    % We also need to check for NaNs; these indicate that an arbitrary
    % coordinate can be chosen
    p(isnan(p)) = 0;
end

% ------------------------------ END OF CODE ------------------------------
