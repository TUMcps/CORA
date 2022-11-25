function [obj]=inputTie(obj,options)
% inputTie - tie: time interval error; computes the error done by 
% the linear assumption of the constant input solution
%
% Syntax:  
%    [obj]=inputTie(obj,options)
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

% Author: Matthias Althoff
% Written:      22-June-2009
% Last update:  13-November-2017
% Last revision:---

%------------- BEGIN CODE --------------

%obtain powers and convert them to interval matrices
Apower=cell(1,length(obj.power.int));
%matrix zonotope
for i=1:length(obj.power.zono)
    Apower{i}=intervalMatrix(obj.power.zono{i});
end
%interval matrix
for i=(length(obj.power.zono)+1):length(obj.power.int)
    Apower{i}=obj.power.int{i};
end

r=obj.stepSize;

%initialize Asum
infimum=-0.25*r;%*r^2;already included in Apower
supremum=0;
timeInterval=intervalMatrix(0.5*(supremum+infimum),0.5*(supremum-infimum));
Asum=timeInterval*Apower{1}*(1/factorial(2));

for i=3:obj.taylorTerms
    %compute factor
    exp1=-i/(i-1); exp2=-1/(i-1);
    infimum = (i^exp1-i^exp2)*r;%*r^i;already included in Apower
    supremum = 0;   
    timeInterval=intervalMatrix(0.5*(supremum+infimum),0.5*(supremum-infimum));
    %compute powers
    Aadd=timeInterval*Apower{i-1};
    %compute sum
    Asum=Asum+Aadd*(1/factorial(i));
end

% compute error due to finite Taylor series according to "M. L. Liou. A novel 
% method of evaluating transient response. In Proceedings of the IEEE, 
% volume 54, pages 20-23, 1966".
% Consider that the power of A is less than for t due to the constant input
% solution instead of the initial state solution
norm_A = norm(obj.A, inf);
epsilon = norm_A*r/(obj.taylorTerms + 2);
if epsilon < 1
    %phi = norm_A^(obj.taylorTerms+1)*r^(obj.taylorTerms+2)/(factorial(obj.taylorTerms + 2)*1-epsilon);
    phi = norm_A^(obj.taylorTerms+1)*r^(obj.taylorTerms+1)/(factorial(obj.taylorTerms + 1)*(1-epsilon));
    Einput = ones(length(obj.A))*interval(-1,1)*phi/norm_A;
else
    Einput = interval(-inf,inf);
    disp('Taylor order not high enough');
end
%write to object structure
obj.inputF=Asum+Einput;
%obj.inputF=Asum+obj.E*r;


%------------- END OF CODE --------------