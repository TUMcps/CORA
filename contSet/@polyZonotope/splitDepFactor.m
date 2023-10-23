function pZsplit = splitDepFactor(pZ,ind,varargin)
% splitDepFactor - Splits one dependent factor of a polynomial zonotope
%
% Syntax:
%    pZsplit = splitDepFactor(pZ,ind)
%    pZsplit = splitDepRactor(pZ,ind,polyOrd)
%
% Inputs:
%    pZ - polyZonotope object
%    ind - identifier of the dependent factor that is splitted
%    polyOrd - maximum number of polynomial terms that are splitted exactly
%              (without an over-approximation)
%
% Outputs:
%    pZsplit - cell array of split polyZonotopes
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1],[1;2]);
%    pZsplit = splitDepFactor(pZ,1,5);
%   
%    plot(pZ,[1,2],'FaceColor','r');
%
%    figure; hold on;
%    plot(pZsplit{1},[1,2],'FaceColor','b');
%    plot(pZsplit{2},[1,2],'FaceColor','g');
%
% Reference:
%   [1] Kochdumper, Niklas. Extensions of Polynomial Zonotopes and their 
%       Application to Verification of Cyber-Physical Systems. Diss. 
%       Technische Universität München, 2022.
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: split, splitLongestGen

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       24-March-2018
% Last update:   19-July-2022 (TL, optimizations)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% find selected dependent factor
ind = pZ.id == ind; % 
if sum(ind) ~= 1
    throw(CORAerror('CORA:wrongValue','second',...
        "Given value for 'ind' should contained in identifiers of polynomial zonotope"));
end

E = pZ.E;
E_ind = E(ind, :);

% parse input arguments
if nargin == 3
    polyOrd = varargin{1}; 
else
    polyOrd = max(E_ind);
end

% determine all generators in which the selected dependent factor occurs
genInd = 0 < E_ind & E_ind <= polyOrd;

% [1, Prop 3.1.43/44] using bounds [0,1] and [-1,0]
% create coeffs for i=1...polyOrd:
% (0.5 + 0.5 * a_ind)^i and (-0.5 + 0.5 * a_ind)^i
polyCoeff1 = cell(max(2,polyOrd),1);
polyCoeff2 = cell(max(2,polyOrd),1);
hout = sum(~genInd);

% create pascal triangle
P = {[1, 1]};
for i=2:polyOrd
    P{i} = [1 sum(P{i-1}([1:(i-1); 2:i])) 1];
end

for i=1:polyOrd
    Pi = P{i};
    polyCoeff1{i} = 0.5^i * Pi;
    polyCoeff2{i} = 0.5^i * Pi .* (-mod(i:-1:0, 2)*2+1);

    numExpi = sum(E_ind == i);
    hout = hout + length(Pi) * numExpi;
end

% construct the modified generators for the splitted zonotopes
c1 = pZ.c;
c2 = pZ.c;

G1 = nan(length(c1), hout);
G2 = nan(length(c2), hout);
Eout = nan(size(E, 1), hout); % identical for splitted sets

h = 1;
dh = sum(~genInd);
G1(:, h:h+dh-1) = pZ.G(:, ~genInd);
G2(:, h:h+dh-1) = pZ.G(:, ~genInd);
Eout(:, h:h+dh-1) = E(:, ~genInd);
h = h + dh;

for i = 1:polyOrd
    coef1 = polyCoeff1{i};
    coef2 = polyCoeff2{i};

    expi = E_ind == i;

    dh = length(coef1) * sum(expi);
    
    G1(:, h:h+dh-1) = kron(coef1, pZ.G(:, expi));
    G2(:, h:h+dh-1) = kron(coef2, pZ.G(:, expi));

    Eout(:, h:h+dh-1) = repmat(E(:, expi),1,i+1);
    Eout(ind, h:h+dh-1) = kron(0:i, ones(1, sum(expi))); % fix a_ind

    h = h + dh;

    genInd = xor(genInd, expi);
end

% over-approximate all selected generators that did not get splitted
if sum(genInd) > 0
    Eout = [Eout;zeros(1,size(Eout,2))];
    
    Eout(end, genInd) = Eout(ind, genInd);
    Eout(ind, genInd) = 0;
end

% add every generator with all-zero exponent matrix to the zonotope center
temp = sum(Eout,1);
genInd = temp == 0;

c1 = c1 + sum(G1(:,genInd),2);
G1(:,genInd) = [];
Eout(:,genInd) = [];

c2 = c2 + sum(G2(:,genInd),2);
G2(:,genInd) = [];

% construct the resulting polynomial zonotopes
pZsplit = cell(1, 2);
pZsplit{1} = polyZonotope(c1,G1,pZ.GI,Eout);
pZsplit{2} = polyZonotope(c2,G2,pZ.GI,Eout);

% ------------------------------ END OF CODE ------------------------------
