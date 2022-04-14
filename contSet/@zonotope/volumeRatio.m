function ratio = volumeRatio(Z,P,dims)
% volume - computes the approximate volume ratio of a zonotope and its
%    over-approximating polytope
%
% Syntax:  
%    ratio = volumeRatio(varargin)
%
% Inputs:
%    Z - zonotope object
%    P - polytope object
%    dims - considered dimensions for the approximation
%
% Outputs:
%    ratio - approximated normalized volume ratio
%
% Example:
%    Z = zonotope([1;0],rand(2,5));
%    P = enclosingPolytope(Z);
%    ratio = volumeRatio(Z,P,1);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       12-September-2008 
% Last update:   28-Aug-2019
% Last revision: ---

%------------- BEGIN CODE --------------

%write inputs to variables
if nargin == 1
    error("Not enough input values");
elseif nargin == 2
    dims=dim(Z);
end

%obtain dimension
n=dim(Z);
%generate dim vector
dimVector=1:dims;
%obtain number of iterations
iterations=n-dims+1;

%init projected zonotope
Zproj=Z;

partialRatio = zeros(iterations,1);
for i=1:iterations
    %projected dimensions
    projDims=dimVector+i-1;
    %project zonotope
    Zproj.Z=Z.Z(projDims,:);
    %project polytope
    Pproj=projection(P,projDims);
    
    %compute volume of the projected zonotope and polytope
    volZ=volume(Zproj);
    volP=volume(Pproj);
    
    %obtain normalized ratio
    partialRatio(i)=(volP/volZ)^(1/dims);
end

%final ratio is the mean value of the partial ratios
ratio=mean(partialRatio);

%------------- END OF CODE --------------