function error_zono = linError_mixed_noInt(nlnsysDT,Rdelta,params,options)
% linError_mixed_noInt - computes the linearization error
%
% Syntax:
%    error_zono = linError_mixed_noInt(nlnsysDT,Rdelta,params,options)
%
% Inputs:
%    nlnsysDT - nonlinearSysDT object
%    Rdelta - initial set
%    params - model parameters
%    options - options struct
%
% Outputs:
%    error_zono - zonotope overapproximating the linearization error
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] M. Althoff and B. Krogh, "Reachability Analysis of Nonlinear
%       Differential-Algebraic Systems" in IEEE Transactions on Automatic 
%       Control, vol. 59, no. 2, pp. 371-383, 2014.

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       21-August-2012
% Last update:   25-July-2016 (intervalhull replaced by interval)
%                29-January-2018 (NK)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% if nonlinearARX: system dimension = dimension of output
if isa(nlnsysDT,"nonlinearARX")
    n = nlnsysDT.dim_y;
else
    n = nlnsysDT.nrOfStates;
end

% use correct Hessian file
nlnsysDT = setHessian(nlnsysDT,'int');

%obtain intervals and combined interval z
dx = interval(Rdelta);
du = interval(params.U - center(params.U));
dz = [dx; du];

%compute interval of reachable set
totalInt_x = dx + nlnsysDT.linError.p.x;

%compute intervals of input
totalInt_u = du + nlnsysDT.linError.p.u;

%compute zonotope of state and input
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

%obtain hessian tensor
if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
   ~strcmp(options.lagrangeRem.method,'interval')

    % create taylor models or zoo-objects
    [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

    % evaluate the hessian tensor 
    H = nlnsysDT.hessian(objX,objU);
else
    H = nlnsysDT.hessian(totalInt_x, totalInt_u);
end

%obtain absolute values
dz_abs = max(abs(infimum(dz)), abs(supremum(dz)));

%separate evaluation
H_mid = cell(n,1);
H_rad = cell(n,1);
for i=1:n
    H_mid{i} = sparse(center(H{i}));
    H_rad{i} = sparse(rad(H{i}));
end
error_mid = 0.5*quadMap(Z,H_mid);

%interval evaluation
error_rad = zeros(n,1);
for i=1:n
    error_rad(i,1) = 0.5*dz_abs'*H_rad{i}*dz_abs;
end

%combine results
error_rad_zono = zonotope(interval(-error_rad,error_rad));
error_zono = error_mid + error_rad_zono;

if contains(options.alg,'adaptive')
    error_zono = reduce(error_zono,'adaptive',options.redFactor);
else
    error_zono = reduce(error_zono,options.reductionTechnique,options.zonotopeOrder);
end

% ------------------------------ END OF CODE ------------------------------
