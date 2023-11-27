function Y = outputSet(obj,options,R)
% outputSet - calculates output set based on a (non-)linear output equation
%
% Syntax:
%    Y = outputSet(obj,options,R)
%
% Inputs:
%    obj - contDynamics object
%    options - options for the computation of reachable sets
%    R - reachable set (either time point [i] or time interval [i,i+1])
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
%                23-June-2023 (LL, consider inputs in first-order term)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% skip computation of output set
if ~options.compOutputSet
    Y = R; return
end

% linParamSys or linProbSys
if isa(obj,'linParamSys') || isa(obj,'linProbSys')
    % ...adapt this when these classes also support output functions
    Y = R; return
end

% nonlinear systems: check tensor order for evaluation of output equation
if ~any(options.tensorOrderOutput == [2,3])
    throw(CORAerror('CORA:notSupported',...
        'Only tensor orders 2 and 3 supported for computation of output set.'));
end

% only nonlinear systems from here on...

if ~iscell(R)
    Y = aux_outputSet(obj,options,R,obj.out_isLinear);
else
    Y = cell(length(R),1);
    for i=1:length(R)
        if isstruct(R{i})
            % time-point solution
            Y{i}.set = aux_outputSet(obj,options,R{i}.set,obj.out_isLinear);
            Y{i}.prev = R{i}.prev;
            if isfield(R{i},'parent')
                Y{i}.parent = R{i}.parent;
            end
        else
            % time-interval solution
            Y{i} = aux_outputSet(obj,options,R{i},obj.out_isLinear);
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function Y = aux_outputSet(obj,options,R,all_linear)
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
p = [p_x;p_u];
if isa(obj,'nonlinParamSys')
    if isa(options.paramInt,'interval')
        p_p = center(options.paramInt);
    else
        p_p = options.paramInt;
    end
end

% interval enclosure of set
I_x = interval(R);
I_u = interval(U);
I = cartProd(I_x,I_u);

% evaluate reset function and Jacobian at expansion point
D_lin = zeros(obj.nrOfOutputs,obj.nrOfInputs);
if isa(obj,'nonlinearSys') || isa(obj,'nonlinearSysDT')
    zerothorder = obj.out_mFile(p_x,p_u);
    [J,D_lin] = obj.out_jacobian(p_x,p_u);
elseif isa(obj,'nonlinParamSys')
    if isa(options.paramInt,'interval')
        % ...copied from nonlinParamSys/linearize
        % constant
        f0cell = obj.out_parametricDynamicFile(p_x,p_u);
        % normalize cells
        for i=2:length(f0cell)
            f0cell{1} = f0cell{1} + center(options.paramInt(i-1))*f0cell{i};
            f0cell{i} = rad(options.paramInt(i-1))*f0cell{i};
        end
        % create constant input zonotope
        f0Mat(length(f0cell{1}(:,1)),length(f0cell)) = 0;
        for i=1:length(f0cell)
            f0Mat(:,i) = f0cell{i};
        end
        zerothorder = zonotope(f0Mat);

        % Jacobian
        J = obj.out_jacobian(p_x,p_u);
        % create matrix zonotopes, convert to interval matrix for Jacobian
        J = intervalMatrix(matZonotope(J{1},J(2:end)));

    else % options.paramInt is numeric
        zerothorder = obj.out_mFile(p_x,p_u,p_p);
        J = obj.out_jacobian_freeParam(p_x,p_u,p_p);
    end
end

% first-order
firstorder = J * (R + (-p_x)) + D_lin * (U + (-p_u));

if all_linear
    % only affine map
    secondorder = 0;
    thirdorder = 0;

elseif options.tensorOrderOutput == 2

    % assign correct hessian (using interval arithmetic)
    obj = setOutHessian(obj,'int');

    % evaluate Hessian using interval arithmetic
    if isa(obj,'nonlinParamSys')
        H = obj.out_hessian(I_x,I_u,options.paramInt);
    else
        H = obj.out_hessian(I_x,I_u);
    end

    % obtain maximum absolute values within Z
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
    if isa(obj,'nonlinParamSys')
        H = obj.out_hessian(p_x,p_u,p_p);
    else
        H = obj.out_hessian(p_x,p_u);
    end

    % Cartesian product of given set of states and inputs
    Z = cartProd(R,I_u);

    % second-order
    secondorder = 0.5*quadMap(Z+(-p),H);
    
    % evaluate third-order tensors over entire set
    [T,ind] = obj.out_thirdOrderTensor(I_x,I_u);
    
    % compute Lagrange remainder
    I = I - p;
    thirdorder = interval(zeros(r,1));
    
    for i = 1:r
        for j = 1:n+m
            if ~isempty(T{i,j}) % same as: ~isempty(ind(i,j)) ?
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
