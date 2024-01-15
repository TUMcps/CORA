function p = randPoint_(hyp,N,type,varargin)
% randPoint_ - generates a random point within a constrained hyperplane
%
% Syntax:
%    p = randPoint_(hyp)
%    p = randPoint_(hyp,N)
%    p = randPoint_(hyp,N,type)
%
% Inputs:
%    hyp - conHyperplane object
%    N - number of random points
%
% Outputs:
%    p - random point or matrix of random points
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],2));
%    p = randPoint(hyp);
%
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(hyp,[1,2]);
%    scatter(p(1,:),p(2,:),16,'r','filled');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/randPoint, interval/randPoint

% Authors:       Victor Gassmann
% Written:       10-June-2022
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename randPoint_)

% ------------------------------ BEGIN CODE -------------------------------

% 'type' standard supported
if ~strcmp(type,'standard')
    throw(CORAerror('CORA:notSupported',...
        "Only type = 'standard' supported for class conHyperplane."));
end

% parse input arguments
if ~(isempty(hyp.C) && all(hyp.d==0))
    throw(CORAerror('CORA:notSupported',...
        'Only implemented for hyperplanes without constraints.'));
end

% get properties (a*x=b)
a = hyp.a';
b = hyp.b;
n = length(a);

% normalization
v = a/norm(a);
c_ = b/norm(a);

% get v'*(x-c) = 0: v'*x = c_ <=> v'(x-v*c_) = 0
c = v*c_;

% compute null space
M = null(v');

% number of vectors in N is equal to n-1
if size(M,2) ~= n-1
    throw(CORAerror('CORA:specialError','Number of null space vectors wrong!'));
end

% initialize matrix of points
Dx = zeros(n,N);
for i=1:N
    % get random linear combination
    facs = randn(1,n-1);
    Dx(:,i) = sum(facs.*M,2);
end

% add center to each sampled point
p = Dx + c;

% ------------------------------ END OF CODE ------------------------------
