function [Y, Verror] = outputSet(sys,R,params,options)
% outputSet - computes output set based on a (non-)linear output equation
%
% Syntax:
%    Y = outputSet(sys,R,params,options)
%
% Inputs:
%    sys - contDynamics object
%    R - reachable set (either time point [i] or time interval [i,i+1])
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    Y - output set (either time point [i] or time interval [i,i+1])
%    Verror - linearization error
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
%                23-June-2023 (LL, consider inputs in first-order term)
%                06-November-2023 (LL, add Verror as output)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% skip computation of output set
if ~options.compOutputSet
    Y = R; return
end

% linParamSys or linProbSys
if isa(sys,'linParamSys') || isa(sys,'linProbSys')
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
    [Y, Verror] = aux_outputSet(sys,options.tensorOrderOutput,R,params);
else
    Y = cell(length(R),1);
    for i=1:length(R)
        if isstruct(R{i})
            % time-point solution
            [Y{i}.set, Verror] = aux_outputSet(sys,...
                options.tensorOrderOutput,R{i}.set,params);
            Y{i}.prev = R{i}.prev;
            if isfield(R{i},'parent')
                Y{i}.parent = R{i}.parent;
            end
        else
            % time-interval solution
            [Y{i}, Verror] = aux_outputSet(sys,...
                options.tensorOrderOutput,R{i},params);
        end
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [Y, Verror] = aux_outputSet(sys,tensorOrderOutput,R,params)
% output set computation for a given reachable set

% dimension of state and inputs
n = sys.nrOfDims;
m = sys.nrOfInputs;
r = sys.nrOfOutputs;

% input set
U = params.U;

% expansion points
p_x = center(R);
p_u = center(U);
p = [p_x;p_u];
if isa(sys,'nonlinParamSys')
    if isa(params.paramInt,'interval')
        p_p = center(params.paramInt);
    else
        p_p = params.paramInt;
    end
end

% interval enclosure of set
I_x = interval(R);
I_u = interval(U);
I = cartProd(I_x,I_u);

% evaluate reset function and Jacobian at expansion point
D_lin = zeros(sys.nrOfOutputs,sys.nrOfInputs);
if isa(sys,'nonlinearSys') || isa(sys,'nonlinearSysDT')
    zerothorder = sys.out_mFile(p_x,p_u);
    [J,D_lin] = sys.out_jacobian(p_x,p_u);
elseif isa(sys,'nonlinParamSys')
    if isa(params.paramInt,'interval')
        % ...copied from nonlinParamSys/linearize
        % constant
        f0_3D = sys.out_parametricDynamicFile(p_x,p_u);
        f0cell = arrayfun(@(i) f0_3D(:,:,i),1:sys.nrOfParam+1,...
            'UniformOutput',false);
        % normalize cells
        for i=2:length(f0cell)
            f0cell{1} = f0cell{1} + center(params.paramInt(i-1))*f0cell{i};
            f0cell{i} = rad(params.paramInt(i-1))*f0cell{i};
        end
        % create constant input zonotope
        f0Mat(length(f0cell{1}(:,1)),length(f0cell)) = 0;
        for i=1:length(f0cell)
            f0Mat(:,i) = f0cell{i};
        end
        zerothorder = zonotope(f0Mat);

        % Jacobian
        J = sys.out_jacobian(p_x,p_u);
        % create matrix zonotopes, convert to interval matrix for Jacobian
        J = intervalMatrix(matZonotope(J(:,:,1),J(:,:,2:end)));

    else % params.paramInt is numeric
        zerothorder = sys.out_mFile(p_x,p_u,p_p);
        J = sys.out_jacobian_freeParam(p_x,p_u,p_p);
    end
end

% first-order
firstorder = J * (R + (-p_x)) + D_lin * (U + (-p_u));

if sys.out_isLinear
    % only affine map
    secondorder = zonotope(zeros(size(J,1),1));
    thirdorder = 0;

elseif tensorOrderOutput == 2

    % assign correct hessian (using interval arithmetic)
    sys = setOutHessian(sys,'int');

    % evaluate Hessian using interval arithmetic
    if isa(sys,'nonlinParamSys')
        H = sys.out_hessian(I_x,I_u,params.paramInt);
    else
        H = sys.out_hessian(I_x,I_u);
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
    thirdorder = zeros(r,1);

elseif tensorOrderOutput == 3

    % set handles to correct files
    sys = setOutHessian(sys,'standard');
    sys = setOutThirdOrderTensor(sys,'int');

    % evaluate Hessians at expansion point
    if isa(sys,'nonlinParamSys')
        H = sys.out_hessian(p_x,p_u,p_p);
    else
        H = sys.out_hessian(p_x,p_u);
    end

    % Cartesian product of given set of states and inputs
    Z = cartProd(R,I_u);

    % second-order
    secondorder = 0.5*quadMap(Z+(-p),H);
    
    % evaluate third-order tensors over entire set
    [T,ind] = sys.out_thirdOrderTensor(I_x,I_u);
    
    % compute Lagrange remainder
    I = I - p;
    thirdorder = interval(zeros(r,1));
    
    for i = 1:r
        for j = 1:n+m
            if ~representsa(T{i,j},'emptySet') % same as: ~isempty(ind(i,j)) ?
                thirdorder(i) = thirdorder(i) ...
                    + I(j) * transpose(I) * T{i,j} * I;
            end
        end
    end
    
    % include factor and convert remainder to zonotope for quadMap below
    thirdorder = zonotope(1/6*thirdorder);

end

% compute output set
Verror = secondorder + thirdorder;
Y = zerothorder + firstorder + Verror;

end

% ------------------------------ END OF CODE ------------------------------
