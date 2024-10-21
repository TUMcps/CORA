function Verr = linError_thirdOrder(nlnsysDT,Rdelta,params,options)
% linError_thirdOrder - computes the linearization error
%
% Syntax:
%    Verr = linError_thirdOrder(nlnsysDT,options,R)
%
% Inputs:
%    nlnsysDT - nonlinearSysDT system object
%    R - actual reachable set
%    params - model parameters
%    options - options struct
%
% Outputs:
%    Verr - set of abstraction errors 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       21-August-2012
% Last update:   25-July-2016 (intervalhull replaced by interval)
%                29-January-2018 (NK)
%                08-April-2021 (NK, removed separated eval. of Lag. rem)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set correct tensor files
nlnsysDT = setHessian(nlnsysDT,'standard');
nlnsysDT = setThirdOrderTensor(nlnsysDT,'int');

% compute interval enclosure of reachable set
dx = interval(Rdelta);
du = interval(params.U - center(params.U));
dz = [dx;du];

Int_x = dx + nlnsysDT.linError.p.x;
Int_u = du + nlnsysDT.linError.p.u;

% reduce order before quadMap to save computation time
if contains(options.alg,'adaptive')
    Rred = reduce(Rdelta,'adaptive',options.redFactor);
else
    Rred = reduce(Rdelta,options.reductionTechnique,options.errorOrder);
end
if isa(nlnsysDT, 'nonlinearARX') && (isa(Rred,'polyZonotope') || ...
            isa(Rred,'conPolyZono'))
    Z = stack(Rred,params.U - center(params.U));
else
    Z = cartProd(Rred,params.U - center(params.U));
end

% calculate hessian tensor
H = nlnsysDT.hessian(nlnsysDT.linError.p.x,nlnsysDT.linError.p.u);

% evaluate third-order tensor with range bounding
if isfield(options,'lagrangeRem') && ...
   isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(Int_x,Int_u,options);

    % evaluate third order tensor 
    [T,ind] = nlnsysDT.thirdOrderTensor(objX,objU);
else
    [T,ind] = nlnsysDT.thirdOrderTensor(Int_x,Int_u);
end

% second order abstraction error
error_secondOrder = 0.5*quadMap(Z,H);

% Lagrange remainder
rem = interval(zeros(nlnsysDT.nrOfStates,1),zeros(nlnsysDT.nrOfStates,1));
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
