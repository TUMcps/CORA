function [Zcubic] = cubicMultiplication(Z,C)
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
% Last update:  11-August-2011
% Last revision:---

%------------- BEGIN CODE --------------

%retrieve dimension and number of generators
dim = length(Z.Z(:,1));
gens = length(Z.Z(1,:))-1;
eqs = length(C(:,1));

%get generators including the center
G = Z.Z;

%compute center
for i=1:eqs
    res_tmp = 0;
    for j=1:dim
        for k=1:dim
            for l=1:dim
                res_other_tmp = 0;
                for s=1:gens
                    res_other_tmp = res_other_tmp + 0.5*C{i,j}(k,l)*(...
                        G(j,s+1)*G(k,s+1)*G(l,1) + ...
                        G(j,s+1)*G(k,1)*G(l,s+1) + ...
                        G(j,1)*G(k,s+1)*G(l,s+1));
                end
                res_tmp = res_tmp + res_other_tmp + ...
                    C{i,j}(k,l)*G(j,1)*G(k,1)*G(l,1);
            end
        end
    end
    c(i,1) = res_tmp;
end


%compute 1st set of generators
G1new = zeros(eqs,gens);
for i=1:eqs
    for s=1:gens
        G1new(i,s) = 0;
        for j=1:dim
            G1new(i,s) = G1new(i,s) + G(:,s+1)'*C{i,j}*G(:,s+1)*G(j,s+1);
        end
    end
end


%compute 2nd set of generators
G2new = zeros(eqs,gens);
for i=1:eqs
    for s=1:gens
        G2new(i,s) = 0;
        for j=1:dim
            G2new(i,s) = G2new(i,s) + 0.5*(...
                G(:,s+1)'*C{i,j}*G(:,1)*G(j,s+1) + ...
                G(:,1)'*C{i,j}*G(:,s+1)*G(j,s+1) + ...
                G(:,s+1)'*C{i,j}*G(:,s+1)*G(j,1));
        end
    end
end


%compute 3rd set of generators
for i=1:length(C(:,1))
    counter = 1;
    for s=0:(gens-1)
        for t=(s+1):gens
            G3new(i,counter) = 0;
            for j=1:dim
                G3new(i,counter) = G3new(i,counter) + ...
                    G(:,s+1)'*C{i,j}*G(:,t+1)*G(j,s+1) + ...
                    G(:,t+1)'*C{i,j}*G(:,s+1)*G(j,s+1) + ...
                    G(:,s+1)'*C{i,j}*G(:,s+1)*G(j,t+1);
            end
            counter = counter + 1;
        end
    end
end

%compute 4th set of generators
for i=1:length(C(:,1))
    counter = 1;
    for s=1:(gens-1)
        for t=(s+1):gens
            G4new(i,counter) = 0;
            for j=1:dim
                G4new(i,counter) = G4new(i,counter) + ...
                    G(:,t+1)'*C{i,j}*G(:,s+1)*G(j,t+1) + ...
                    G(:,s+1)'*C{i,j}*G(:,t+1)*G(j,t+1) + ...
                    G(:,t+1)'*C{i,j}*G(:,t+1)*G(j,s+1);
            end
            counter = counter + 1;
        end
    end
end


%compute 5th set of generators
for i=1:length(C(:,1))
    counter = 1;
    for s=0:(gens-2)
        for t=(s+1):(gens-1)
            for u=(t+1):gens
                G5new(i,counter) = 0;
                for j=1:dim
                    G5new(i,counter) = G5new(i,counter) + ...
                        G(:,t+1)'*C{i,j}*G(:,u+1)*G(j,s+1) + ...
                        G(:,u+1)'*C{i,j}*G(:,t+1)*G(j,s+1) + ...
                        G(:,s+1)'*C{i,j}*G(:,u+1)*G(j,t+1) + ...
                        G(:,s+1)'*C{i,j}*G(:,u+1)*G(j,t+1) + ...
                        G(:,u+1)'*C{i,j}*G(:,s+1)*G(j,t+1) + ...
                        G(:,s+1)'*C{i,j}*G(:,t+1)*G(j,u+1); 
                      
                end
                counter = counter + 1;
            end
        end
    end
end


%generate new zonotope
Zcubic = zonotope([c, G1new, G2new, G3new, G4new, G5new]);

%------------- END OF CODE --------------