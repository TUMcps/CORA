function [obj]=tie(obj)
% tie - tie: time interval error; computes the error done by 
% building the convex hull of time point solutions
%
% Syntax:  
%    [obj]=tie(obj)
%
% Inputs:
%    obj - linear interval system object
%    options - options struct
%
% Outputs:
%    obj - linear interval system object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: expm, inputSol

% Author:       Matthias Althoff
% Written:      06-August-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain powers and convert them to interval matrices
Apower=cell(1,length(obj.power.int));
%matrix zonotope
for i=2:length(obj.power.zono)
    Apower{i}=intervalMatrix(obj.power.zono{i});
end
%interval matrix
for i=(length(obj.power.zono)+1):length(obj.power.int)
    Apower{i}=obj.power.int{i};
end

%initialize Asum
inf=-0.25;%*r^2;already included in Apower
sup=0;
timeInterval=intervalMatrix(0.5*(sup+inf),0.5*(sup-inf));
Asum=timeInterval*Apower{2}*(1/factorial(2));

for i=3:obj.taylorTerms
    %compute factor
    exp1=-i/(i-1); exp2=-1/(i-1);
    inf = (i^exp1-i^exp2);%*r^i;already included in Apower
    sup = 0;   
    timeInterval=intervalMatrix(0.5*(sup+inf),0.5*(sup-inf));
    %compute powers
    Aadd=timeInterval*Apower{i};
    %compute sum
    Asum=Asum+Aadd*(1/factorial(i));
end

%write to object structure
obj.F=Asum+obj.E;

%------------- END OF CODE --------------