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
%    P = sumPoints(randn(2,4),randn(2,5));
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

if isempty(varargin)
    throw(CORAerror('CORA:notEnoughInputArgs',1));
end

if ~ismatrix(varargin{1}) 
    throw(CORAerror('CORA:wrongValue','first','Must be double matrix'));
end
n = size(varargin{1},1);

for i=1:length(varargin)
    inputArgsCheck({{varargin{i},'att',{'double'},{'nrows',n}}});
end

V = combineVec(varargin{:});
N = length(varargin);
n = size(varargin{1},1);
P = zeros(n,size(V,2));
for i=1:n
    P(i,:) = sum(V(i:n:N*n,:),1);
end
P = unique(P','stable','rows')';

% ------------------------------ END OF CODE ------------------------------
