function Y = outputSet(obj,options,R,R_y)
% outputSet - calculates output set based on a (non-)linear output equation
%
% Syntax:
%    Y = outputSet(obj,options,R,R_y)
%
% Inputs:
%    obj - nonlinDASys object
%    options - options for the computation of reachable sets
%    R - reachable set (either time point [i] or time interval [i,i+1])
%    R_y - reachable set of algebraic variables (either time point [i] or
%          time interval [i,i+1])
%
% Outputs:
%    Y - output set (either time point [i] or time interval [i,i+1])
%
% Example:
%    -
%
% References:
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

if ~iscell(R)
    Y = aux_outputSet(obj,options,R,R_y,obj.out_isLinear);
else
    Y = cell(length(R),1);
    for i=1:length(R)
        if isstruct(R{i})
            % time-point solution
            Y{i}.set = aux_outputSet(obj,options,R{i}.set,R_y{i},obj.out_isLinear);
            Y{i}.prev = R{i}.prev;
        else
            % time-interval solution
            Y{i} = aux_outputSet(obj,options,R{i},R_y{i},obj.out_isLinear);
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function Y = aux_outputSet(obj,options,R,R_y,all_linear)
% output set computation for a given reachable set

% dimension of state and inputs
n = obj.dim;
m = obj.nrOfInputs;
r = obj.nrOfOutputs;

% input set
U = options.U + options.uTrans;

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
zerothorder = obj.out_mFile(p_x,p_y,p_u);

% evaluate Jacobian at expansion point
J = obj.out_jacobian(p_x,p_y,p_u);

% first-order
firstorder = J * (R + (-p_x));

if all_linear
    % only affine map
    secondorder = 0;
    thirdorder = 0;

elseif options.tensorOrderOutput == 2

    % assign correct hessian (using interval arithmetic)
    obj = setOutHessian(obj,'int');

    % evaluate Hessian using interval arithmetic
    H = obj.out_hessian(I_x,I_y,I_u);

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

elseif options.tensorOrderOutput == 3

    % set handles to correct files
    obj = setOutHessian(obj,'standard');
    obj = setOutThirdOrderTensor(obj,'int');

    % evaluate Hessians at expansion point
    H = obj.out_hessian(p_x,p_y,p_u);

    % Cartesian product of given set of states and inputs
    Z = cartProd(R,[I_y;I_u]);

    % second-order
    secondorder = 0.5*quadMap(Z + (-p),H);
    
    % evaluate third-order tensors over entire set
    T = obj.out_thirdOrderTensor(I_x,I_y,I_u);
    
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
