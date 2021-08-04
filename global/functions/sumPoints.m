function P = sumPoints(varargin)
% sumPoints - Computes all possible sums of given input arguments
%
% Syntax:  
%    P = sumPoints(varargin)
%
% Inputs:
%    ~ - any number of double matrices
%
% Outputs:
%    P - sum of points 
%
% Example: 
%    P = sumPoints(randn(2,4),randn(4,5));
%
% Other m-files required: -
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      22-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if isempty(varargin)
    error('No input arguments provided!');
end

% check if all are doubles
if ~all(cellfun(@(inp)isa(inp,'double'),varargin))
    error('Some inputs are not of type double!');
end

% check if all inputs have same row-dimension
if length(unique(cellfun(@(inp) size(inp,1),varargin)))~=1
    error('Not all inputs have same row dimension!');
end

V = combineVec(varargin{:});
N = length(varargin);
n = size(varargin{1},1);
P = zeros(n,size(V,2));
for i=1:n
    P(i,:) = sum(V(i:n:N*n,:),1);
end
P = unique(P','stable','rows')';
%------------- END OF CODE --------------