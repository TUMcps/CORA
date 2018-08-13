function obj = intersectionTimeModel(obj,options)
% intersectionTimeModel - creates a polynomial model to determine the time 
% when a linear system hits an arbitrary switching surface h(x)
%
% Syntax:  
%    [timeModelFile, jacobianFile] = intersectionTimeModel(obj)
%
% Inputs:
%    obj - switchingSurface object
%
% Outputs:
%
%    timeModelFile - handle to the time model file
%    jacobianFile - handle to the jacobian of the time model file
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      21-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%obtain symmetric varibles for the intersection time model
[x_0,x_c,A,u,t_0,a_s,B_s,C_s,y] = symVariables(obj, 'LRbrackets');

%set number of parameters and store y
obj.nrOfParameters = length(y);
obj.y_sym = y;

%create model
if obj.timeOrder == 1
    t_s = t_0 + a_s.'*(x_0 - x_c);
elseif obj.timeOrder == 2
    t_s = t_0 + a_s.'*(x_0 - x_c) + (x_0 - x_c).'*B_s*(x_0 - x_c);
elseif obj.timeOrder == 3
    %todo
    t_s = [];
end

%create trajectory (input u not yet considered)
truncated_eAt = eye(obj.dim);
for iOrder = 1:obj.stateOrder
    truncated_eAt = truncated_eAt + (A*t_s)^iOrder/factorial(iOrder);
end
x = truncated_eAt*x_0;

%insert into switching surface equation
eq = obj.mFile(x);

%write eq into mFile
createIntersectionTimeModelFile(eq,options.path,func2str(obj.mFile));

%compute jacobian with respect to various parameters
J = jacobian(eq,y);

%write jacobian to file
createJacobianFile(J,options.path,func2str(obj.mFile));





%------------- END OF CODE --------------