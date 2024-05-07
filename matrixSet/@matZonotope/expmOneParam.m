function [eZ,eI,zPow,iPow,E,RconstInput] = expmOneParam(matZ,r,maxOrder,varargin)
% expmOneParam - operator for the exponential matrix of a matrix zonotope,
%    evaluated dependently. Higher order terms are computed via interval
%    arithmetic.
%
% Syntax:
%    [eZ,eI,zPow,iPow,E,RconstInput] = expmOneParam(matZ,r,maxOrder,varargin)
%
% Inputs:
%    matZ - matZonotope object
%    r - time step size
%    intermediate Order - Taylor series order until computation is 
%                         performed with matrix zonotopes
%    maxOrder - maximum Taylor series order until remainder is computed
%    options - options struct
%
% Outputs:
%    eZ - matrix zonotope exponential part
%    eI - interval matrix exponential part
%    zPow - ?
%    iPow - cell array storing the powers of the matrix:
%           A,A^2,...,A^(intermediateOrder)
%    E - interval matrix for the remainder
%    RconstInput - ???
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       13-September-2010 
% Last update:   04-April-2017
%                26-April-2024 (TL, using new matZonotope class)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%cannot directly use u as input since zonotope has preference over
%matZonotopes
if length(varargin) == 1
    options = varargin{1};
    if ~isa(options.uTrans,'zonotope')
        u = zonotope([options.uTrans,zeros(size(options.uTrans))]);
    else
        u = zonotope(options.uTrans);
    end
else
    u = zonotope([0,0]);
end

%obatin matrix center and generator
C = matZ.center;
G1 = matZ.G(:,:,1); % TL, correct? Was matZ.generators{1} before

%obtain center and generator of input uTrans
c_u = u.c;
g_u = u.G;

% pre-allocate matrices D,E (center and generators of powers)
D = nan([size(C),maxOrder]);
E = nan([size(G1),maxOrder,maxOrder]);
if isscalar(c_u)
    % use dimension of state
    D_u = nan(size(D));
    E_u = nan(size(E));
else
    % use dimension of input
    D_u = nan([size(c_u),maxOrder]);
    E_u = nan([size(g_u),maxOrder,maxOrder]);
end
zPow = cell(1,maxOrder);

% init first entry
D(:,:,1) = C;
E(:,:,1,1) = G1;

D_u(:,:,1) = c_u;
E_u(:,:,1,1) = g_u;

%update power cell
zPow{1} = matZ*r; 

%the first cell index refers to the power!
for n = 2:maxOrder
    D(:,:,n) = D(:,:,n-1)*C;
    E(:,:,1,n) = D(:,:,n-1)*G1 + E(:,:,1,n-1)*C;
    for i = 2:(n-1)
        E(:,:,i,n) = E(:,:,i-1,n-1)*G1 + E(:,:,i,n-1)*C;
    end
    E(:,:,n,n) = E(:,:,n-1,n-1)*G1;
    
    %input
    D_u(:,:,n) = D(:,:,n-1)*c_u;
    E_u(:,:,1,n) = D(:,:,n-1)*g_u + E(:,:,1,n-1)*c_u;
    for i = 2:(n-1)
        E_u(:,:,i,n) = E(:,:,i-1,n-1)*g_u + E(:,:,i,n-1)*c_u;
    end
    E_u(:,:,n,n) = E(:,:,n-1,n-1)*g_u;
    
    %update power cell
    zPow{n} = matZonotope(D(:,:,n),E(:,:,1:n,n))*r^n; 
end

% compute exponential matrix

% preallocate
E_sum = nan([size(E,1:2),maxOrder]);
E_u_sum = nan([size(E_u,1:2),maxOrder]);

%generators
for n = 1:maxOrder
    factor = r^n/factorial(n);
    E_sum(:,:,n) = E(:,:,1,n)*factor;
    E_u_sum(:,:,n) = E_u(:,:,1,n)*factor;
    for i=(n+1):maxOrder
        factor = r^i/factorial(i);
        E_sum(:,:,n) = E_sum(:,:,n) + E(:,:,n,i)*factor;
        E_u_sum(:,:,n) = E_u_sum(:,:,n) + E_u(:,:,n,i)*factor;
    end
end

%center
D_sum = eye(dim(matZ)) + D(:,:,1)*r;
D_u_sum = D_u(:,:,1)*r;
for i = 2:maxOrder
    factor = r^i/factorial(i);
    D_sum = D_sum + D(:,:,i)*factor;
    D_u_sum = D_u_sum + D_u(:,:,i)*factor;
end

%reduce size of generators for even numbers and add to center
for n = 1:floor(maxOrder/2)
    E_sum(:,:,2*n) = 0.5*E_sum(:,:,2*n);
    D_sum = D_sum + E_sum(:,:,2*n);
    
    E_u_sum(:,:,2*n) = 0.5*E_u_sum(:,:,2*n);
    D_u_sum = D_u_sum + E_u_sum(:,:,2*n);
end

%compute remainder
matI = intervalMatrix(matZ*r);
E = exponentialRemainder(matI,maxOrder);

%write result to eZ and eI
eZ = matZonotope(D_sum, E_sum);
eI = E;

%obtain constant input zonotope
RconstInput = zonotope(matZonotope(D_u_sum, E_u_sum));
%RconstInput = zonotope(D_u_sum, E_u_sum);

% no powers based on interval matrix
iPow = [];

% ------------------------------ END OF CODE ------------------------------
