function Verr = linError_thirdOrder(obj, options, R)
% linError_thirdOrder - computes the linearization error
%
% Syntax:
%    Verr = linError_thirdOrder(obj,options,R)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    Verr - set of abstraction errors 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       21-August-2012
% Last update:   25-July-2016 (intervalhull replaced by interval)
%                29-January-2018 (NK)
%                08-April-2021 (NK, removed separated eval. of Lag. rem)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set correct tensor files
obj = setHessian(obj,'standard');
obj = setThirdOrderTensor(obj,'int');

% compute interval enclosure of reachable set
dx = interval(R);
du = interval(options.U);
dz = [dx;du];

Int_x = dx + obj.linError.p.x;
Int_u = du + obj.linError.p.u;

% reduce order before quadMap to save computation time
if contains(options.alg,'adaptive')
    Rred = reduce(R,'adaptive',options.redFactor);
else
    Rred = reduce(R,options.reductionTechnique,options.errorOrder);
end
Z = cartProd(Rred,options.U);

% calculate hessian tensor
H = obj.hessian(obj.linError.p.x,obj.linError.p.u);

% evaluate third-order tensor with range bounding
if isfield(options,'lagrangeRem') && ...
   isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(Int_x,Int_u,options);

    % evaluate third order tensor 
    [T,ind] = obj.thirdOrderTensor(objX,objU);
else
    [T,ind] = obj.thirdOrderTensor(Int_x,Int_u);
end

% second order abstraction error
error_secondOrder = 0.5*quadMap(Z,H);

% Lagrange remainder
rem = interval(zeros(obj.dim,1),zeros(obj.dim,1));
for i=1:length(ind)
    temp = interval(0,0);
    for j=1:length(ind{i})
        temp = temp + (dz.' * T{i,ind{i}(j)}*dz) * dz(ind{i}(j));
    end
    rem(i,1) = 1/6*temp;
end

% overall abstraction error
Verr = error_secondOrder + rem;
if contains(options.alg,'adaptive')
    Verr = reduce(Verr,'adaptive',options.redFactor);    
else
    Verr = reduce(Verr,options.reductionTechnique,options.zonotopeOrder);    
end

% ------------------------------ END OF CODE ------------------------------
