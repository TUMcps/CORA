function ratio = volumeRatio(Z,P,dims)
% volumeRatio - computes the approximate volume ratio of a zonotope and its
%    over-approximating polytope
%
% Syntax:
%    ratio = volumeRatio(Z,P,dims)
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
%    P = polytope(Z);
%    ratio = volumeRatio(Z,P,1)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       12-September-2008 
% Last update:   28-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%write inputs to variables
if nargin == 1
    throw(CORAerror('CORA:notEnoughInputArgs',2));
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
    Zproj = project(Z, projDims);
    %project polytope
    Pproj=project(P,projDims);
    
    %compute volume of the projected zonotope and polytope
    volZ=volume(Zproj);
    volP=volume(Pproj);
    
    %obtain normalized ratio
    partialRatio(i)=(volP/volZ)^(1/dims);
end

%final ratio is the mean value of the partial ratios
ratio=mean(partialRatio);

% ------------------------------ END OF CODE ------------------------------
