function [c_int, deltaVec] = cubicMultiplication_interval(Z,Cint)
% cubicMultiplication - computes \{C_{ilmn}*x_l*x_l*x_n|x \in Z\}
%
% Syntax:  
%    [Zcubic] = cubicMultiplication_interval(Z,C)
%
% Inputs:
%    Z - zonotope object
%    C - 4 dimensional tensor
%
% Outputs:
%    Zcubic - zonotope object
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
% Written:      11-August-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%retrieve dimension and number of generators
gens = length(Z.Z(1,:))-1;
dim = length(Z.Z(:,1));

%get center and generators
c = Z.Z(:,1);
G = Z.Z(:,2:end);

%square computation--------------------------------------------------------
%new center
for i = 1:dim
    c_int(i,1) = interval(0,0);
    for l = 1:dim
        for m = 1:dim
            for n = 1:dim
                centerCorr = 0;
                for o = 1:gens
                    centerCorr = centerCorr + 0.5*(c(n)*G(m,o)*G(l,o) + c(l)*G(n,o)*G(m,o) + c(l)*G(m,o)*G(n,o));
                end
                c_int(i,1) = c_int(i,1) + Cint{i}(l,m,n)*(c(l)*c(m)*c(n) + centerCorr);
            end
        end
    end
end

%define delta vector
deltaVec = zeros(dim,1);

%new linear part
iGen = 0;
for o = 1:gens
    iGen = iGen+1;
    for i = 1:dim
        lin_int(i,iGen) = interval(0,0);
        for l = 1:dim
            for m = 1:dim
                for n = 1:dim
                    lin_int(i,iGen) = lin_int(i,iGen) + Cint{i}(l,m,n)...
                        *(c(n)*c(m)*G(l,o) + c(l)*c(n)*G(m,o) + c(l)*c(m)*G(n,o));
                end
            end
        end
    end
end

%add to delta vector
for i=1:iGen
    deltaVec = deltaVec + supremum(abs(lin_int(:,iGen)));
end

%new quad-linear part
iGen = 0;
for o = 1:gens
    iGen = iGen+1;
    for i = 1:dim
        quad_lin_int(i,iGen) = interval(0,0);
        for l = 1:dim
            for m = 1:dim
                for n = 1:dim
                    quad_lin_int(i,iGen) = quad_lin_int(i,iGen) + Cint{i}(l,m,n)...
                        *(c(n)*G(m,o)*G(l,o) + c(l)*G(n,o)*G(m,o) + c(l)*G(m,o)*G(n,o));
                end
            end
        end
    end
end

%add to delta vector
for i=1:iGen
    deltaVec = deltaVec + supremum(abs(quad_lin_int(:,iGen)));
end

%new quad part
iGen = 0;
for o = 1:(gens-1)
    for p = (o+1):gens
        iGen = iGen+1;
        for i = 1:dim
            quad_int(i,iGen) = interval(0,0);
            for l = 1:dim
                for m = 1:dim
                    for n = 1:dim
                        quad_int(i,iGen) = quad_int(i,iGen) + Cint{i}(l,m,n)...
                            *(c(n)*G(m,o)*G(l,p) + c(l)*G(n,o)*G(m,p) + c(l)*G(m,o)*G(n,p));
                    end
                end
            end
        end
    end
end

%add to delta vector
for i=1:iGen
    deltaVec = deltaVec + supremum(abs(quad_int(:,iGen)));
end


%new cubic part
iGen = 0;
for o = 1:gens
    for p = 1:gens
        for q = 1:gens
            iGen = iGen+1;
            for i = 1:dim
                cubic_int(i,iGen) = interval(0,0);
                for l = 1:dim
                    for m = 1:dim
                        for n = 1:dim
                            cubic_int(i,iGen) = cubic_int(i,iGen) + Cint{i}(l,m,n)*G(l,o)*G(m,p)*G(n,q);
                        end
                    end
                end
            end
        end
    end
end

%add to delta vector
for i=1:iGen
    deltaVec = deltaVec + supremum(abs(cubic_int(:,iGen)));
end


%------------- END OF CODE --------------