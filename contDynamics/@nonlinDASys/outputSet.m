function Y = outputSet(nlnsysDA,R,R_y,params,options)
% outputSet - calculates output set based on a (non-)linear output equation
%
% Syntax:
%    Y = outputSet(nlnsysDA,R,R_y,params,options)
%
% Inputs:
%    nlnsysDA - nonlinDASys object
%    R - reachable set (either time point [i] or time interval [i,i+1])
%    R_y - reachable set of algebraic variables (either time point [i] or
%          time interval [i,i+1])
%    params - model parameters (U, uTrans)
%    options - options for the computation of reachable sets
%
% Outputs:
%    Y - output set (either time point [i] or time interval [i,i+1])
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-November-2022
% Last update:   07-December-2022 (MW, allow to skip output set)
%                30-August-2024 (MW, add params to input arguments)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% skip computation of output set
if ~options.compOutputSet
    Y = R; return
end

% check tensor order
if ~any(options.tensorOrderOutput == [2,3])
    throw(CORAerror('CORA:notSupported',...
        'Only tensor orders 2 and 3 supported for computation of output set.'));
end

% compute input
U = params.U + params.uTrans;

if ~iscell(R)
    Y = aux_outputSet(nlnsysDA,options.tensorOrderOutput,R,R_y,U);
else
    Y = cell(length(R),1);
    for i=1:length(R)
        if isstruct(R{i})
            % time-point solution
            Y{i}.set = aux_outputSet(nlnsysDA,options.tensorOrderOutput,...
                R{i}.set,R_y{i},U);
            Y{i}.prev = R{i}.prev;
        else
            % time-interval solution
            Y{i} = aux_outputSet(nlnsysDA,options.tensorOrderOutput,R{i},R_y{i},U);
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function Y = aux_outputSet(nlnsysDA,tensorOrderOutput,R,R_y,U)
% output set computation for a given reachable set

% dimension of state and inputs
n = nlnsysDA.nrOfDims;
m = nlnsysDA.nrOfInputs;
r = nlnsysDA.nrOfOutputs;

% expansion points
p_x = center(R);
p_u = center(U);
p_y = center(R_y);
p = [p_x;p_y;p_u];

% interval enclosure of set
I_x = interval(R);
I_u = interval(U);
I_y = interval(R_y);
I = [I_x;I_y;I_u];

% evaluate reset function at expansion point
zerothorder = nlnsysDA.out_mFile(p_x,p_y,p_u);

% evaluate Jacobian at expansion point
J = nlnsysDA.out_jacobian(p_x,p_y,p_u);

% first-order
firstorder = J * (R + (-p_x));

if nlnsysDA.out_isLinear
    % only affine map
    secondorder = 0;
    thirdorder = 0;

elseif tensorOrderOutput == 2

    % assign correct hessian (using interval arithmetic)
    nlnsysDA = setOutHessian(nlnsysDA,'int');

    % evaluate Hessian using interval arithmetic
    H = nlnsysDA.out_hessian(I_x,I_y,I_u);

    % obtain maximum absolute values within I
    dz = max(abs(infimum(I)),abs(supremum(I)));

    % calculate the second-order error
    secondorder = zeros(length(H),1);
    for i = 1:length(H)
        H_ = abs(H{i});
        H_ = max(infimum(H_),supremum(H_));
        secondorder(i) = 0.5 * dz' * H_ * dz;
    end
    secondorder = zonotope(zeros(r,1),diag(secondorder));

    % no third-order computation
    thirdorder = 0;

elseif tensorOrderOutput == 3

    % set handles to correct files
    nlnsysDA = setOutHessian(nlnsysDA,'standard');
    nlnsysDA = setOutThirdOrderTensor(nlnsysDA,'int');

    % evaluate Hessians at expansion point
    H = nlnsysDA.out_hessian(p_x,p_y,p_u);

    % Cartesian product of given set of states and inputs
    Z = cartProd(R,[I_y;I_u]);

    % second-order
    secondorder = 0.5*quadMap(Z + (-p),H);
    
    % evaluate third-order tensors over entire set
    T = nlnsysDA.out_thirdOrderTensor(I_x,I_y,I_u);
    
    % compute Lagrange remainder
    I = I - p;
    thirdorder = interval(zeros(r,1));
    
    for i = 1:r
        for j = 1:n+m
            if ~representsa_(T{i,j},'emptySet',eps)
                thirdorder(i) = thirdorder(i) ...
                    + I(j) * transpose(I) * T{i,j} * I;
            end
        end
    end
    
    % include factor and convert remainder to zonotope for quadMap below
    thirdorder = zonotope(1/6*thirdorder);

end

% compute output set
Y = zerothorder + firstorder + secondorder + thirdorder;

end

% ------------------------------ END OF CODE ------------------------------
