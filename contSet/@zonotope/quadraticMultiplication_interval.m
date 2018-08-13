function [Zquad] = quadraticMultiplication_interval(Z,Q)
% quadraticMultiplication - computes \{Q_{ijk}*x_j*x_k|x \in Z\}
%
% Syntax:  
%    [Zquad] = quadraticMultiplication(Z1,Q)
%
% Inputs:
%    Z - zonotope object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    Zquad - zonotope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      18-May-2011
% Last update:  07-July-2017
% Last revision:---

%------------- BEGIN CODE --------------

%retrieve dimension and number of generators
gens = length(Z.Z(1,:))-1;

%get center and generators
c = Z.Z(:,1);
G = Z.Z(:,2:end);

%square computation--------------------------------------------------------
%new center
sq_c = sqMult(Q,c);
for i=1:gens
    quad(:,i) = sqMult(Q,G(:,i));
    sq_c = sq_c + 0.5*quad(:,i);
end

%get generators
sq_G = interval();
%1st set of generators
for i=1:gens
    sq_G(:,end+1)=mixedMult(Q,c,G(:,i)) + mixedMult(Q,G(:,i),c);
end
%2nd set of generators
for i=1:gens
    sq_G(:,end+1)=0.5*quad(:,i);
end
%get indices for 3rd set of generators
if (gens>=2)
    ind = combinator(gens,2,'c');
    for i=1:length(ind(:,1))
        ind1 = ind(i,1);
        ind2 = ind(i,2);
        sq_G(:,end+1) = mixedMult(Q,G(:,ind1),G(:,ind2))...
                     + mixedMult(Q,G(:,ind2),G(:,ind1));
    end
end
%--------------------------------------------------------------------------

%separate uncertain center and generators in certain ones
sq_c_cer = mid(sq_c);
err = rad(sq_c);

for i=1:length(sq_G(1,:))
    sq_G_cer(:,i) = mid(sq_G(:,i));
    err = err + rad(sq_G(:,i));
end

%generate new zonotope
Zquad = zonotope([sq_c_cer, sq_G_cer, diag(err)]);

%delete zeros
Zquad=deleteZeros(Zquad);

function y = sqMult(Q,x)

dim = length(Q);
for i=1:dim
    y(i,1) = interval(x'*Q{i}*x);
end

function y = mixedMult(Q,x1,x2)

dim = length(Q);
for i=1:dim
    y(i,1) = interval(x1'*Q{i}*x2);
end

%------------- END OF CODE --------------