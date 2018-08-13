function X_0 = plotPointCloud(obj,A,u,x_c,R_0,dims)
% plotPointCloud - plot point cloud of solutions according to the
% polynomial intersection time model
%
% Syntax:  
%    plotPointCloud(obj,A,u,R_0)
%
% Inputs:
%    obj - switchingSurface object
%    A - system maztrix
%    u - constant input
%    R_0 - initial set
%
% Outputs:
%    ---
%
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      22-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%retrieve time model parameters
[t_0, a_s, B_s, C_s] = getIntersectionTimeParameters(obj);

%obtain random points from initial set
nrOfPoints = 100;
%X_0 = pointSetExtreme(R_0,nrOfPoints);
X_0 = pointSet(R_0,nrOfPoints);

%compute solution for each point (input not yet considered!)
for i = 1:nrOfPoints
    x_curr = X_0(:,i);
    if obj.timeOrder == 1
        t_sol = t_0 + a_s.'*(x_curr - x_c);
    elseif obj.timeOrder == 2
        t_sol = t_0 + a_s.'*(x_curr - x_c) + (x_curr - x_c).'*B_s*(x_curr - x_c);
    elseif obj.timeOrder == 3
        cubic_sol = [];
        t_sol = t_0 + a_s.'*(x_curr - x_c) + (x_curr - x_c).'*B_s*(x_curr - x_c);
    end
    %exact solution
    x_sol = expm(A*t_sol)*x_curr;
    
    %approximative solution
    expmApprox = eye(obj.dim);
    for i=1:obj.stateOrder
        expmApprox = expmApprox + (A*t_sol)^i/factorial(i);
    end
    x_approx = expmApprox*x_curr;
    
    %plot
    plot(x_sol(dims(1)), x_sol(dims(2)), '.r');
    plot(x_approx(dims(1)), x_approx(dims(2)), 'ok');
end
 

%------------- END OF CODE --------------