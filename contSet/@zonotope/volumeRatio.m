function [ratio] = volumeRatio(varargin)
% volume - Computes the approxiamte volume ratio of a zonotope and its
% overapproximating polytope
%
% Syntax:  
%    [ratio] = volumeRatio(varargin)
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
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 12-September-2008 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%write inputs to variables
if nargin==2
    Z=varargin{1};
    P=varargin{2};
    dims=length(Z.Z(:,1));    
elseif nargin==3
    Z=varargin{1};
    P=varargin{2};
    dims=varargin{3};
end

%obtain dimension
dim=length(Z.Z(:,1));
%generate dim vector
dimVector=1:dims;
%obtain number of iterations
iterations=dim-dims+1;

%init projected zonotope
Zproj=Z;

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