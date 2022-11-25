function [obj] = nonlinMap(obj,options)
% nonlinMap - computes the interval matrix which abstracts the nonlinear
% map due to the higher order terms in the Taylor expansion
%
% Syntax:  
%    [obj] = nonlinMap(obj,options)
%
% Inputs:
%    obj - linear interval system object
%    options - options struct
%
% Outputs:
%    obj - linear interval system object
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      03-January-2009 
% Last update:  22-June-2009
%               22-July-2009
%               03-February-2010
% Last revision: ---

%------------- BEGIN CODE --------------

    
%load data from object/options structure
A=obj.A;
taylorTerms=options.taylorTerms;
r=options.timeStep;
dim=dimension(obj);

%compute exact A square
Asquare=apprSquare(A);
Aint=intervalMatrix(A);

%initialize overapproximation
Apower{1}=Aint;
Asum=0*Apower{1};

%compute powers for each term and sum of these
for i=3:taylorTerms
    %compute powers
    Apower{i-1}=Apower{i-2}*Aint;
    
    %compute sums
    Asum=Asum+Apower{i-2}*r^(i-2)/factorial(i);
end
Asum=Asquare*r^2*(0.5*eye(dim)+Asum);

%add missing powers
for i=(taylorTerms-1):taylorTerms
    %compute powers
    Apower{i+1}=Apower{i}*Aint;    
end

%determine error due to finite Taylor series
alpha=infNorm(Aint);
epsilon=alpha*r/(taylorTerms+2);
phi=(alpha*r)^(taylorTerms+1)/factorial(taylorTerms+1)/(1-epsilon);
E=interval(-ones(dim),ones(dim))*phi;

%nonlinear map
N=Asum+E;

%write to object structure
obj.taylor.N=N; %nonlinear map N
obj.taylor.powers=Apower;
obj.taylor.error=E;    
obj.taylor.sq=Asquare;


%------------- END OF CODE --------------