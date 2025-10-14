function P = sumPoints(varargin)
% sumPoints - Computes all possible sums of given input arguments
%
% Syntax:
%    P = sumPoints(varargin)
%
% Inputs:
%    ~ - any number of double matrices with same number of rows
%
% Outputs:
%    P - sum of points 
%
% Example: 
%    P1 = [ -0.359 0.734 1.136 0.472 ; -0.130 0.120 -0.687 0.288 ];
%    P2 = [ 1.392 0.001 -2.342 2.809 0.287 ; -1.346 0.053 1.248 -0.232 -0.465 ];
%    S = sumPoints(P1,P2);
%
% Other m-files required: -
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       22-March-2021
% Last update:   17-November-2022 (VG, reworked inputArgsCheck)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,Inf);
if ~ismatrix(varargin{1}) 
    throw(CORAerror('CORA:wrongValue','first','Must be double matrix'));
end
n = size(varargin{1},1);
for i=1:length(varargin)
    inputArgsCheck({{varargin{i},'att','double',@(value) size(value,1) == n}});
end

% compute
V = combineVec(varargin{:});
N = length(varargin);
n = size(varargin{1},1);
P = zeros(n,size(V,2));
for i=1:n
    P(i,:) = sum(V(i:n:N*n,:),1);
end
P = unique(P','stable','rows')';

% ------------------------------ END OF CODE ------------------------------
