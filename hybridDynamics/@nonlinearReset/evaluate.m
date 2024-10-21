function x_ = evaluate(nonlinReset,x,varargin)
% evaluate - evaluates the nonlinear reset function for a given state (set)
%    and input (set)
%
% Syntax:
%    x_ = evaluate(nonlinReset,x)
%    x_ = evaluate(nonlinReset,x,u)
%
% Inputs:
%    nonlinReset - nonlinearReset object
%    x - state before reset
%    u - (optional) input before reset
%
% Outputs:
%    x_ - state after reset
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linReset/evaluate

% Authors:       Mark Wetzlinger
% Written:       13-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default value
narginchk(2,3);
u = setDefaultValues({zeros(nonlinReset.inputDim,1)},varargin);
inputArgsCheck({{nonlinReset,'att','nonlinearReset','scalar'};...
                {x,'att',{'contSet','numeric'}};...
                {u,'att',{'contSet','numeric'}}});

% special case: no sets (can be extend to representsa_(...,'point',eps) but
% potentially time-consuming...?)
if isnumeric(x) && isnumeric(u)
    x_ = nonlinReset.f(x,u);
    return
end

% check if derivatives have been computed...
if isempty(nonlinReset.J)
    throw(CORAerror('CORA:specialError','Derivatives have not been computed.'));
    % filepath = [CORAROOT filesep 'models' filesep 'auxiliary' filesep 'nonlinearReset' filesep];
    % fname = 'reset_function';
    % nonlinReset = derivatives(nonlinReset,filepath,fname,2);
end

% expansion point for Taylor series
[px,pu] = aux_expansionPoint(x,u);
z = aux_cat_state_input(x,u);

% interval conversion/enclosure of sets
Ix = interval(x);
Iu = interval(u);

% evaluate reset function at expansion point
f = nonlinReset.f(px,pu);

% evaluate Jacobian at expansion point
[Jx,Ju] = nonlinReset.J(px,pu);

if nonlinReset.tensorOrder == 2
    [secondOrder,thirdOrder] = aux_tensorOrder2(nonlinReset,px,pu,Ix,Iu);
else  %tensorOrder == 3
    throw(CORAerror('CORA:notSupported',...
        'nonlinearReset/evaluate only supports tensor order 2.'));
    [secondOrder,thirdOrder] = aux_tensorOrder3(nonlinReset,z,px,pu,Ix,Iu);
end

% evaluate enclosure
x_ = f + [Jx,Ju]*(z-[px;pu]) + secondOrder + thirdOrder;

end


% Auxiliary functions -----------------------------------------------------

function [px,pu] = aux_expansionPoint(x,u)
% expansion point is chosen at the center of the set

    if isa(x,'contSet')
        px = center(x);
    elseif isnumeric(x)
        px = x;
    end
    
    if isa(u,'contSet')
        pu = center(u);
    elseif isnumeric(x)
        pu = u;
    end

end

function z = aux_cat_state_input(x,u)

    if isa(x,'contSet') || isa(u,'contSet')
        z = cartProd(x,u);
    else  % both are numeric
        z = [x;u];
    end

end

function [secondOrder,thirdOrder] = aux_tensorOrder2(nonlinReset,px,pu,Ix,Iu)

% dimensions
n = nonlinReset.preStateDim;
m = nonlinReset.inputDim;

% evaluate Hessian over entire set
[Hx,Hu] = nonlinReset.H(Ix,Iu);
% concatenate by hand... (i always 1:2 because of *two* variables x, u)
H = arrayfun(@(i) [interval(reshape(Hx.inf(1,:,:),[n+m,n]), ...
                            reshape(Hx.sup(1,:,:),[n+m,n])), ...
                   interval(reshape(Hu.inf(1,:,:),[n+m,m]), ...
                            reshape(Hu.sup(1,:,:),[n+m,m]))], ...
             1:2,'UniformOutput',false);

% obtain maximum absolute values within sets for x and u shifted by
% their respective linearization points
dx = max(abs(infimum(Ix - px)),abs(supremum(Ix - px)));
du = max(abs(infimum(Iu - pu)),abs(supremum(Iu - pu)));
dz = [dx;du];

secondOrder = zeros(nonlinReset.postStateDim,1);
for i=1:numel(H)
    H_abs = abs(H{i});
    H_ = max(infimum(H_abs),supremum(H_abs));
    secondOrder(i) = 0.5 * dz' * H_ * dz;
end
secondOrder = interval(-secondOrder,secondOrder);

thirdOrder = 0;

end

function [secondOrder,thirdOrder] = aux_tensorOrder3(nonlinReset,z,px,pu,Ix,Iu)

% evaluate Hessians at expansion point
[Hx,Hu] = nonlinReset.H(px,pu);
% rewrite as cell for quadratic map
H = aux_reshapeHessian(Hx,Hu,n,m);
% Q = repmat({0},[n,1]);
% for i = 1:n
%     funHan = trans.reset.Q{i};
%     if ~isnumeric(funHan)
%         % if Hessian is numeric, it is zero
%         Q{i} = funHan(p);
%     end
% end
secondOrder = 0.5*quadMap(z,H);

% evaluate third-order tensors over entire set
[Tx,Tu] = nonlinReset.T(Ix,Iu);
% T = repmat({0},[n,n+m]);
% for i = 1:n
%     for j = 1:n+m
%         funHan = trans.reset.T{i,j};
%         if ~isnumeric(funHan)
%             % if third-order tensor is numeric, it is zero
%             T{i,j} = funHan(I);
%         end
%     end
% end
Ix_ = Ix - px;
Iu_ = Iu - pu;
L = interval(zeros(nonlinReset.postStateDim,1));
for i = 1:n
    for j = 1:n+m
        % L(i) = L(i) + Ix_(j) * transpose(Ix_) * Tx{i,j} * Ix_;
        % L(i) = L(i) + Iu_(j) * transpose(Iu_) * Tu{i,j} * Iu_;
    end
end

% include factors from derivatives
thirdOrder = 1/6*L;

end

function H = aux_reshapeHessian(Hx,Hu,n,m)
% concatenate Hessian tensors so that they can be used in quadMap operation

% allocate cell array
H = cell(n+m,1);

% convert state Hessians
if n > 1
    for i=1:n
        H{i} = interval(Hx.inf(:,:,i),Hx.sup(:,:,i));
    end
else
    H{1} = Hx;
end

% convert input Hessians
if m > 1
    for i=1:m
        H{n+i} = interval(Hu.inf(:,:,i),Hu.sup(:,:,i));
    end
else
    H{end} = Hu;
end

end

% ------------------------------ END OF CODE ------------------------------
