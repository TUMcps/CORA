function [P]=randomPolytope(dim,points,type,maxLength)
% randomPolytope - generates a random polytope 
%
% Syntax:  
%    [P]=randomPolytope(dim,points,type)
%
% Inputs:
%    dim - dimension
%    points - number of points
%    type - selector for different probability distributions of the
%    point to center distance
%
% Outputs:
%    P - polytope object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      01-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%generate random vector for the generator lengths
%uniform distribution
if strcmp(type,'uniform')
    for i=1:points
        %uniform distribution
        l(i) = unifrnd(0,maxLength);
    end
    %     l=rand(generators,1);
    
%exponential distribution    
elseif strcmp(type,'exponential')
    for i=1:points
        %uniform distribution
        l(i) = exprnd(1);
    end    
%     a=1;
%     tmp=rand(generators,1);
%     l=-log(tmp)/a;

%gaussian distribution    
elseif strcmp(type,'gamma')
    for i=1:points
        %uniform distribution
        l(i) = gamrnd(2,1);
    end    
%     a=1;
%     tmp=rand(generators,1);
%     l=-log(tmp)/a;
end

%create generators
for i=1:points
    %generate random point on sphere
    gTmp=randomPointOnSphere(dim);
    %stretch
    V(i,:)=l(i)*gTmp;
end

%center is set to 0
%create MPT polytope
P = mptPolytope(V);



%------------- END OF CODE --------------