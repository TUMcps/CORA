function [H,Zdelta,errorStat,T,ind3,Zdelta3] = precompStatError_adaptive(obj,options,Rdelta)
% same as precompStatError, but using adaptive zonotope orders

% precompute the second order static error along with hessian matrix
obj = setHessian(obj,'standard');

% initialize output arguments
T = []; ind3 = []; Zdelta3 = [];

% reduce the reachable set for the initial time point
Rred = reduce(Rdelta,'adaptive',sqrt(options.redFactor));

% over-approximation of the reachable set at the initial time point
Rdelta = reduce(zonotope(Rdelta),'adaptive',sqrt(options.redFactor));

% extend the sets by the input sets
Z = cartProd(Rred,options.U);
Zdelta = cartProd(Rdelta,options.U);

% calculate the hessian tensor
H = obj.hessian(obj.linError.p.x, obj.linError.p.u);

% calculate the quadratic map == static second-order error
errorSecOrdStat = 0.5*quadMap(Z,H);

% reduce the order of the set of static errors
errorStat = reduce(errorSecOrdStat,'adaptive',sqrt(options.redFactor));

end