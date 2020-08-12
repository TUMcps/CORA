function [obj] = exponential(obj,options)
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
% Last update:  10-August-2010
%               03-September-2013
% Last revision:---

%------------- BEGIN CODE --------------
   
%load data from object/options structure
A=obj.A;
A_abs=abs(A);
taylorTerms=options.taylorTerms;
r=options.timeStep;
dim=dimension(obj);


%initialize 
Apower{1}=A;  
Apower_abs{1}=A_abs; 
M = eye(dim);
    
%compute powers for each term and sum of these
for i=1:taylorTerms
    %compute powers
    Apower{i+1}=Apower{i}*A;
    Apower_abs{i+1}=Apower_abs{i}*A_abs;
    M = M + Apower_abs{i}*r^(i)/factorial(i);
end   
%determine error due to finite Taylor series
W=expm(A_abs*r)-M;
%compute absolute value of W for numerical stability
W=abs(W);
E=interval(-W,W);
    
%write to object structure
obj.taylor.powers=Apower;
obj.taylor.error=E;    

%------------- END OF CODE --------------