function [min,max] = quadraticProgramming_interval(Z,Q)
% quadraticProgramming - computes \{Q_{ijk}*x_j*x_k|x \in Z\}
%
% Syntax:  
%    [Zquad] = quadraticProgramming_interval(Z,Q)
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
% Written:      31-May-2011
% Last update:  ---
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

%get generators
sq_G = intval();
%1st set of generators
for i=1:gens
    sq_G(:,end+1)=mixedMult(Q,c,G(:,i)) + mixedMult(Q,G(:,i),c);
end
%2nd set of generators
ind = combinator(gens,2,'p','r');
for iDim=1:length(Q)
    for i=1:length(ind(:,1))
        ind1 = ind(i,1);
        ind2 = ind(i,2);
        B{iDim}(ind1,ind2) = G(:,ind1)'*Q{iDim}*G(:,ind2);
    end
end
%--------------------------------------------------------------------------

%random evaluation
[min,max] = randomEval(B,1e4);
min = 0.5*min;
max = 0.5*max;


function y = sqMult(Q,x)

dim = length(Q);
for i=1:dim
    y(i,1) = x'*Q{i}*x;
end

function y = mixedMult(Q,x1,x2)

dim = length(Q);
for i=1:dim
    y(i,1) = x1'*Q{i}*x2;
end


function [min,max] = randomEval(B,samples)

%get dimension
dim = length(B{1});

%initialize
min = zeros(length(B),1);
max = zeros(length(B),1);

for iSample = 1:samples

    %get random vector
    alpha = 1-2*rand(dim,1);

    %evaluation
    for i=1:length(B)
        val = alpha'*B{i}*alpha;
        if sup(val)>max(i)
            max(i) = sup(val);
        elseif inf(val)<min(i)
            min(i) = inf(val);
        end
    end
end



%------------- END OF CODE --------------