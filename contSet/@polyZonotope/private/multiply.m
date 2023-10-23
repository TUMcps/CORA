function [G,E] = multiply(G1,E1,G2,E2)
% multiply - ???
%
% Syntax:
%    [G,E] = multiply(G1,E1,G2,E2)
%
% Inputs:
%    G1 - ???
%    E1 - ???
%    G2 - ???
%    E2 - ???
%
% Outputs:
%    G - ???
%    E - ???
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Victor Gassmann
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if size(G1,2)~=size(E1,2) || size(G2,2)~=size(E2,2)
    if size(G1,2)~=size(E1,2)
        throw(CORAerror('CORA:dimensionMismatch',G1,E1));
    end
    if size(G2,2)~=size(E2,2)
        throw(CORAerror('CORA:dimensionMismatch',G2,E2));
    end
end
if any([size(G1,1),size(G2,1)]~=1)
    throw(CORAerror('CORA:dimensionMismatch',G1,G2));
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
    G = zeros(size(G,1),1);
    E = zeros(size(E,1),1);
end

% ------------------------------ END OF CODE ------------------------------
