function [Z] = parallelotope(varargin)
% parallelotope - Takes the directions in W to enclose a polytope without
% computing the vertices
%
% Syntax:  
%    [Z] = parallelotope(P,W)
%
% Inputs:
%    P - ppl Polytope object
%
% Outputs:
%    Z - parallelotope in zonotope representation
%
% Example: 
%
% Other m-files required: ---
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      26-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%one input
if nargin == 1
    obj = varargin{1};
    W = eye(obj.dim);
%two inputs    
elseif nargin == 2
    obj = varargin{1};
    W = varargin{2};
end

%transform matrix
obj = pinv(W)*obj;

%compute bound matrix
boundMat = bound(-obj.C,obj.d,eye(obj.dim));
lowerLimit = boundMat(:,1);
upperLimit = boundMat(:,2);

%generate center and generators
c=zeros(obj.dim,1);
for i=1:length(lowerLimit)
    %center
    c(i)=0.5*(lowerLimit(i)+upperLimit(i));
    %generator
    Gdiag(i)=0.5*(lowerLimit(i)-upperLimit(i));
end

%generate zonotope
Z=W*zonotope([c,diag(Gdiag)]);





%------------- END OF CODE --------------