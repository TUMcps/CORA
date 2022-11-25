function [G,E] = multiply(G1,E1,G2,E2)
if size(G1,2)~=size(E1,2) || size(G2,2)~=size(E2,2)
    error('sizes of E/G not consistent');
end
if any([size(G1,1),size(G2,1)]~=1)
    error('why?');
end
m1 = size(E1,2);
m2 = size(E2,2);
Et = repmat(E1,1,m2) + repelem(E2,1,m1);
Gt = repmat(G1,1,m2).* repelem(G2,1,m1);
[E,G] = removeRedundantExponents(Et,Gt);
% remove exponents with zero G entries
if isa(G,'double')
    ind = all(G==0,1);
    G(:,ind) = [];
    E(:,ind) = [];
end
% prevent empty result
if isempty(G)
    G = 0;
    E = 0;
end


