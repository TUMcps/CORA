function p = randPoint(hyp,varargin)
% randPoint - generates a random point within a constrained hyperplane
%
% Syntax:  
%    p = randPoint(hyp)
%    p = randPoint(hyp,N)
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
% See also: interval/randPoint

% Author:        Victor Gassmann
% Written:       10-June-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_randPoint('conHyperplane',hyp,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of vars
    p = vars{1}; return
else
    % assign variables
    C = vars{1}; N = vars{2}; type = vars{3}; pr = vars{4};
    % sampling with 'gaussian' is done in contSet method
    if strcmp(type,'gaussian')
        p = randPoint@contSet(C,N,type,pr); return
    end
end

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

% get properties (a'*x=b)
a = hyp.h.c;
b = hyp.h.d;
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

%------------- END OF CODE --------------