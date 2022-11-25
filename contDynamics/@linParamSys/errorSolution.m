function inputSet = errorSolution(obj,options,V)
% errorSolution - computes the solution due to the linearization error
%
% Syntax:  
%    inputSet = errorSolution(obj,options,V)
%
% Inputs:
%    obj - linear system object
%    options - options struct
%    V - set of linearization errors
%
% Outputs:
%    inputSet - reachable set due to the linearization error
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      18-January-2008 
% Last update:  27-May-2011
%               25-July-2016 (intervalhull replaced by interval)
%               01-November-2017 (deletion of zero generators added)
%               04-May-2020 (reordering of input arguments)
% Last revision:---

%------------- BEGIN CODE --------------

%load data from object/options structure
r=options.timeStep;

%initialize the reachable set due to input
inputSet = V*r;

%matrix zonotope
for i=1:length(obj.power.zono)
    inputSet = inputSet + obj.power.zono_input{i}*V;
end
%interval matrix
for i=(length(obj.power.zono)+1):length(obj.power.int)
    taylorTerm = obj.power.int{i}*(r/factorial(i+1));
    inputSet = inputSet + taylorTerm*V;
end

%remainder term
Vabs=zonotope(interval(-1,1)*supremum(abs(interval(V))));
inputSet = inputSet + obj.E*r*Vabs;

%delete zero generators in zonotope representation
inputSet = deleteZeros(inputSet);


%------------- END OF CODE --------------