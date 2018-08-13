function [obj,L,Q,C,L_tilde,Q_tilde] = computeTaylorTensors(obj,A,u,x_c)
% computeTaylorTensors - compute tensors required to map a set of initial
% states onto the switching surface
%
% Syntax:  
%    obj = computeTaylorTensors(obj)
%
% Inputs:
%    obj - switchingSurface object
%
% Outputs:
%    obj - switchingSurface object
%
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      21-August-2013
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%retrieve time model parameters
[t_0, a_s, B_s, C_s] = getIntersectionTimeParameters(obj);

%create auxiliary matrices
matTmp_1 = zeros(obj.dim);
for i=1:obj.stateOrder
    matTmp_1 = matTmp_1 + A^i/factorial(i-1)*t_0^(i-1);
end
matTmp_2 = zeros(obj.dim);
for i=2:obj.stateOrder
    matTmp_2 = matTmp_2 + A^i/(2*factorial(i-2))*t_0^(i-2);
end



%linear part
L_tilde = eye(obj.dim);
for i=1:obj.stateOrder
    L_tilde = L_tilde + A^i/factorial(i)*t_0^i;
end

%quadratic part
%Q tensor
for i = 1:obj.dim
    Q_tilde{i} = a_s*matTmp_1(i,:);
end

%cubic part
%C tensor
for i = 1:obj.dim
    for j = 1:obj.dim
        C_tilde{i,j} = B_s*matTmp_1(i,j) + a_s*a_s.'*matTmp_2(i,j);
    end 
end

%corrections
for i = 1:obj.dim
    for j = 1:obj.dim
        L(i,j) = L_tilde(i,j) - x_c.'*Q_tilde{i}(:,j) + x_c.'*C_tilde{i,j}*x_c;
        Q{i}(j,:) = Q_tilde{i}(j,:) - 2*x_c.'*C_tilde{i,j};
    end
end
C = C_tilde;


%write to object
obj.taylorTensors.L = L;
obj.taylorTensors.Q = Q;
obj.taylorTensors.C = C;

%------------- END OF CODE --------------