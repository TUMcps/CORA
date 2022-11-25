function [eZ,eI,zPow,iPow,E,RconstInput] = expmOneParam(matZ,r,maxOrder,varargin)
% expmMixed - operator for the exponential matrix of a 
% matrix zonotope, evaluated dependently. Higher order terms are computed
% via interval arithmetic.
%
% Syntax:  
%    [eZ,eI,zPow,iPow,E,RconstInput] = expmOneParam(matZ,r,maxOrder,varargin)
%
% Inputs:
%    matZ - matrix zonotope
%    r - time increment
%    intermediate Order - Taylor series order until computation is 
%                           performed with matrix zonotopes
%    maxOrder - maximum Taylor series order until remainder is computed
%    options - options struct
%
% Outputs:
%    eZ - matrix zonotope exponential part
%    eI - interval matrix exponential part
%    zPow - #nodef
%    iPow - cell array storing the powers of the matrix:
%           A,A^2,...,A^(intermediateOrder)
%    E - interval matrix for the remainder
%    RconstInput - #nodef
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      13-September-2010 
% Last update:  04-April-2017
% Last revision:---

%------------- BEGIN CODE --------------

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
G = matZ.generator{1};

%obtain center and generator of input uTrans
u_mat = u.Z;
c_u = u_mat(:,1);
g_u = u_mat(:,2);

%initialize matrices D,E (center and generators of powers)
D{1} = C;
E{1}{1} = G;

D_u{1} = c_u;
E_u{1}{1} = g_u;

%update power cell
zPow{1} = matZ*r; 

%the first cell index refers to the power!
for n = 2:maxOrder
    D{n} = D{n-1}*C;
    E{n}{1} = D{n-1}*G + E{n-1}{1}*C;
    for i = 2:(n-1)
        E{n}{i} = E{n-1}{i-1}*G + E{n-1}{i}*C;
    end
    E{n}{n} = E{n-1}{n-1}*G;
    
    %input
    D_u{n} = D{n-1}*c_u;
    E_u{n}{1} = D{n-1}*g_u + E{n-1}{1}*c_u;
    for i = 2:(n-1)
        E_u{n}{i} = E{n-1}{i-1}*g_u + E{n-1}{i}*c_u;
    end
    E_u{n}{n} = E{n-1}{n-1}*g_u;
    
    %update power cell
    zPow{n} = matZonotope(D{n},E{n})*r^n; 
end

%compute exponential matrix
%generators
for n = 1:maxOrder
    factor = r^n/factorial(n);
    E_sum{n} =  E{n}{1}*factor;
    E_u_sum{n} = E_u{n}{1}*factor;
    for i=(n+1):maxOrder
        factor = r^i/factorial(i);
        E_sum{n} = E_sum{n} + E{i}{n}*factor;
        E_u_sum{n} = E_u_sum{n} + E_u{i}{n}*factor;
    end
end

%center
D_sum = eye(matZ.dim) + D{1}*r;
D_u_sum = D_u{1}*r;
for i = 2:maxOrder
    factor = r^i/factorial(i);
    D_sum = D_sum + D{i}*factor;
    D_u_sum = D_u_sum + D_u{i}*factor;
end

%reduce size of generators for even numbers and add to center
for n = 1:floor(maxOrder/2)
    E_sum{2*n} = 0.5*E_sum{2*n};
    D_sum = D_sum + E_sum{2*n};
    
    E_u_sum{2*n} = 0.5*E_u_sum{2*n};
    D_u_sum = D_u_sum + E_u_sum{2*n};
end

%compute remainder
matI = intervalMatrix(matZ*r);
E = exponentialRemainder( matI,maxOrder);

%write result to eZ and eI
eZ = matZonotope(D_sum, E_sum);
eI = E;

%obtain constant input zonotope
RconstInput = zonotope(matZonotope(D_u_sum, E_u_sum));
%RconstInput = zonotope(D_u_sum, E_u_sum);

   
iPow = []; %no powers based on interval matrix


%------------- END OF CODE --------------