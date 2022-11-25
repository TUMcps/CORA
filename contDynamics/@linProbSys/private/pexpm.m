function [obj] = pexpm(obj,options)
% expm - computes the overapproximation of the exponential of a system 
% matrix up to a certain accuracy
%
% Syntax:  
%    [obj] = expm(obj)
%
% Inputs:
%    obj - linearSys object
%    options - reachability options
%
% Outputs:
%    obj - linearSys object
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      07-May-2007 
% Last update:  08-September-2009
% Last revision:---

%------------- BEGIN CODE --------------
   
%load data from object/options structure
A=obj.A;
taylorTerms=options.taylorTerms;
r=options.timeStep;
dim=dimension(obj);


%initialize 
Apower{1}=A;  
    
%compute powers for each term and sum of these
for i=1:taylorTerms
    %compute powers
    Apower{i+1}=Apower{i}*A;
end   
%determine error due to finite Taylor series
alpha=norm(A,inf);
epsilon=alpha*r/(taylorTerms+2);
phi=(alpha*r)^(taylorTerms+1)/factorial(taylorTerms+1)/(1-epsilon);  
E=interval(-ones(dim),ones(dim))*phi;
    
%write to object structure
obj.taylor.eAt=expm(A*r);
obj.taylor.powers=Apower;
obj.taylor.error=E;      

%------------- END OF CODE --------------