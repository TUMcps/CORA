function [Z_error, errorInt, errorInt_x, errorInt_y, R_y] = linError(obj, options, R, Verror_y)
% linError - computes the linearization error
%
% Syntax:
%    [obj] = linError(obj,options)
%
% Inputs:
%    obj - nonlinear differential algebraic system object
%    options - options struct
%    R - actual reachable set
%    Verror_y - set of algebraic linearization error
%
% Outputs:
%    Z_error - zonotope overapproximating the linearization error
%    errorInt - interval overapproximating the linearization error
%    errorInt_x - interval overapproximating the linearization error (dynamic part)
%    errorInt_y - interval overapproximating the linearization error (constraint part)
%    R_y - reachable set of the algebraic part
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       29-October-2007 
% Last update:   22-January-2008
%                02-February-2010
%                13-February-2012    
%                15-June-2016
%                25-July-2016 (intervalhull replaced by interval)
%                04-August-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set handle to correct file
obj = setHessian(obj,'int');

%compute set of algebraic variables
f0_con = obj.linError.f0_con;
D = obj.linError.D;
E = obj.linError.E;
F_inv = obj.linError.F_inv;
R_y = -F_inv*(f0_con + D*R + E*options.U + Verror_y);

%obtain intervals and combined interval z
dx = interval(R);
dy = interval(R_y);
du = interval(options.U);
dz = [dx; dy; du];

%compute interval of reachable set
totalInt_x = dx + obj.linError.p.x;

%compute interval of algebraic states
totalInt_y = dy + obj.linError.p.y;

%compute intervals of input
totalInt_u = du + obj.linError.p.u;

%obtain Hessian
[Hf, Hg] = obj.hessian(totalInt_x, totalInt_y, totalInt_u);

%error_x
for i=1:length(Hf)
    error_x(i,1) = 0.5*dz.'*Hf{i}*dz;
end

%error_y
for i=1:length(Hg)
    error_y(i,1) = 0.5*dz.'*Hg{i}*dz;
end

% %obtain intervals and combined interval z
% dx = sup(abs(interval(IH_x)));
% dy = sup(abs(interval(IH_y)));
% du = sup(abs(interval(IH_u)));
% dz = [dx; dy; du];
% 
% tic
% %obtain hessian tensor
% [Hf, Hg] = hessianTensor(totalInt_x, totalInt_y, totalInt_u);
% 
% %error_x
% for i=1:length(Hf)
%     error_x(i) = 0.5*dz'*sup(abs(Hf{i}))*dz;
% end
% 
% %error_y
% for i=1:length(Hg)
%     error_y(i) = 0.5*dz'*sup(abs(Hg{i}))*dz;
% end
% toc
% 
% %compute linearization error by passing the intervals to the Lagrange
% %remainder mFile
% tic
% [error_x, error_y] = lagrangeRemainder(totalInt_x,totalInt_y,totalInt_u,dx,dy,du);
% toc

%compute final error
Z_error_x = zonotope(interval(infimum(error_x),supremum(error_x)));
Z_error_y = zonotope(interval(infimum(error_y),supremum(error_y)));
Z_error = Z_error_x + obj.linError.CF_inv*Z_error_y;

%update R_y
R_y =  obj.linError.p.y + (-F_inv)*(f0_con + D*R + E*options.U + Z_error_y);

%error intervals
errorIHabs = abs(interval(Z_error));
errorInt = supremum(errorIHabs);

errorIHabs_y = abs(interval(Z_error_y));
errorInt_y = supremum(errorIHabs_y);

errorIHabs_x = abs(interval(Z_error_x));
errorInt_x = supremum(errorIHabs_x);

% ------------------------------ END OF CODE ------------------------------
