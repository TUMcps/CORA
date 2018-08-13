function [Z]=randomZonotope(dim,generators,type,maxLength)
% randomZonotope - generates a random zonotope 
%
% Syntax:  
%    [Z]=randomZonotope(dim,generators,type)
%
% Inputs:
%    dim - dimension
%    generators - number of generators
%    type - selector for different probability distributions of the
%    generator length
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      30-September-2008
% Last update:  14-February-2011
% Last revision:---

%------------- BEGIN CODE --------------

%generate random vector for the generator lengths
%uniform distribution
if strcmp(type,'uniform')
    for i=1:generators
        %uniform distribution
        %l(i) = unifrnd(0,maxLength);
        l(i) = maxLength*rand(1);
    end
    %     l=rand(generators,1);
    
%exponential distribution    
elseif strcmp(type,'exponential')
    for i=1:generators
        %uniform distribution
        l(i) = exprnd(1);
    end    
%     a=1;
%     tmp=rand(generators,1);
%     l=-log(tmp)/a;

%gaussian distribution    
elseif strcmp(type,'gamma')
    for i=1:generators
        %uniform distribution
        l(i) = gamrnd(2,1);
    end    
%     a=1;
%     tmp=rand(generators,1);
%     l=-log(tmp)/a;
end

%create generators
for i=1:generators
    %generate random point on sphere
    gTmp=randomPointOnSphere(dim);
    %stretch
    G(:,i)=l(i)*gTmp;
end

%center is set 0
c=zeros(dim,1);

%instantiate zonotope
Z=zonotope([c,G]);



%------------- END OF CODE --------------
