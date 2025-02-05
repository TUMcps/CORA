function obj = priv_pexpm(obj,options)
% priv_pexpm - computes the overapproximation of the exponential of a system 
%    matrix up to a certain accuracy
%
% Syntax:
%    obj = priv_pexpm(obj,options)
%
% Inputs:
%    obj - linProbSys object
%    options - reachability options
%
% Outputs:
%    obj - linProbSys object
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       07-May-2007 
% Last update:   08-September-2009
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
   
%load data from object/options structure
A=obj.A;
taylorTerms=options.taylorTerms;
r=options.timeStep;
n=obj.nrOfDims;

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
E=interval(-ones(n),ones(n))*phi;
    
%write to object structure
obj.taylor.eAt=expm(A*r);
obj.taylor.powers=Apower;
obj.taylor.error=E;      

% ------------------------------ END OF CODE ------------------------------
