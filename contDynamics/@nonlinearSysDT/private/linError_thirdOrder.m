function errorZon = linError_thirdOrder(obj, options, R)
% linError_thirdOrder - computes the linearization error
%
% Syntax:  
%    error = linError_thirdOrder(obj,options,R)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    errorZon - zonotope overapproximating the linearization error
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  25-July-2016 (intervalhull replaced by interval)
%               29-January-2018 (NK)
% Last revision:---

%------------- BEGIN CODE --------------

%obtain intervals and combined interval z
dx = interval(R);
du = interval(options.U);
dz = [dx; du];

%compute intervalhull of reachable set
totalInt_x = dx + obj.linError.p.x;

%compute intervals of input
totalInt_u = du + obj.linError.p.u;

%compute zonotope of state and input
Rred = reduce(R,options.reductionTechnique,options.errorOrder);
Z=cartProd(Rred,options.U);

%obtain absolute values
dz_abs = max(abs(infimum(dz)), abs(supremum(dz)));

% calculate hessian tensor
H = obj.hessian(obj.linError.p.x, obj.linError.p.u);

% evaluate third-order tensor
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

    % evaluate third order tensor 
    [T, ind] = obj.thirdOrderTensor(objX, objU);

else
    [T, ind] = obj.thirdOrderTensor(totalInt_x, totalInt_u);
end

%second order error
error_secondOrder = 0.5*quadMap(Z, H);

%interval evaluation
% for i=1:length(T(:,1))
%     error_sum = interval(0,0);
%     for j=1:length(T(1,:))
%         error_tmp = dz'*T{i,j}*dz;
%         error_sum = error_sum + error_tmp * dz(j);
%     end
%     error_thirdOrder_old(i,1) = 1/6*error_sum;
% end
% 
% error_thirdOrder_old_zono = zonotope(error_thirdOrder_old);

%alternative: separate evaluation without interval arithmetic
for i=1:length(T(:,1))
    for j=1:length(T(1,:))
        T_mid{i,j} = sparse(center(T{i,j}));
        T_rad{i,j} = sparse(rad(T{i,j}));
    end
end

error_mid = 1/6*cubMap(Z, T_mid);

%interval evaluation
if isempty(cell2mat(ind))
    % empty zonotope if all entries in T are empty
    error_rad_zono = zonotope(zeros(obj.dim,1));
else
    error_rad = zeros(obj.dim,1);
    for i=1:length(ind)
        error_sum2 = 0;
        for j=1:length(ind{i})
            error_tmp2 = dz_abs'*T_rad{i,j}*dz_abs;
            error_sum2 = error_sum2 + error_tmp2 * dz_abs(j);
        end
        error_rad(i,1) = 1/6*error_sum2;
    end
    error_rad_zono = zonotope(interval(-error_rad, error_rad));
end

%combine results
error_thirdOrder = error_mid + error_rad_zono;
errorZon = error_secondOrder + error_thirdOrder;

errorZon = reduce(errorZon,options.reductionTechnique,options.zonotopeOrder);


%------------- END OF CODE --------------