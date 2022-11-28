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
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: split, splitLongestGen

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      24-March-2018
% Last update:  19-July-2022 (TL: optimizations)
% Last revision:---

%------------- BEGIN CODE --------------

% find selected dependent factor
index = find(pZ.id == ind);
if isempty(index)
    throw(CORAerror('CORA:wrongValue','second',...
        "Given value for 'ind' should contained in identifiers of polynomial zonotope"));
end

% determine all generators in which the selected dependent factor occurs
ind = find(pZ.expMat(index,:) > 0);

% parse input arguments
if nargin == 3
    polyOrd = varargin{1}; 
else
    polyOrd = max(pZ.expMat(index,ind));
end

% create pascal triangle
A = {[1, 1]};
for i=2:polyOrd
    A{i} = [1 sum(A{i-1}([1:(i-1); 2:i])) 1];
end

% create coeffs for (0.5 + 0.5x)^p and (-0.5 + 0.5x)^p
polyCoeff1 = cell(max(2,polyOrd),1);
polyCoeff2 = cell(max(2,polyOrd),1);

for i=1:polyOrd
    Ai = A{i};
    polyCoeff1{i} = 0.5^i * Ai;
    polyCoeff2{i} = 0.5^i * Ai .* (-mod(i:-1:0, 2)*2+1);
end

% construct the modified generators for the splitted zonotopes
expMat1 = pZ.expMat;
% expMat2 = pZ.expMat; % equal to expMat1!
G1 = pZ.G;
G2 = pZ.G;
c1 = pZ.c;
c2 = pZ.c;

for i = 1:polyOrd
    coef1 = polyCoeff1{i};
    coef2 = polyCoeff2{i};

    indTemp = find(pZ.expMat(index,ind) == i);
    n = length(indTemp);
    adv = length(coef1)-1;

    G1_ = zeros(size(G1, 1), n*adv);
    G2_ = zeros(size(G2, 1), n*adv);
    E_ = zeros(size(expMat1, 1), n*adv);

    idx = 1;

    for j = 1:length(indTemp)
       
        g = pZ.G(:,ind(indTemp(j)));
        e = pZ.expMat(:,ind(indTemp(j)));
       
        % first splitted zonotope
        G1(:,ind(indTemp(j))) = coef1(1) * g;
        G1_(:, idx:idx+adv-1) = g*coef1(2:end);

        % expMat is equal for both
        expMat1(index,ind(indTemp(j))) = 0;
        Ej_ = e* ones(1,length(coef1)-1);
        Ej_(index,:) = 1:length(coef1)-1;
        E_(:, idx:idx+adv-1) = Ej_;
        
        % second splitted zonotope
        G2(:,ind(indTemp(j))) = coef2(1) * g;
        G2_(:, idx:idx+adv-1) = g*coef2(2:end);

        idx = idx + adv;
    end

    G1 = [G1, G1_];
    G2 = [G2, G2_];
    expMat1 = [expMat1, E_];
    
    % remove the finished indices from the list
    ind(indTemp) = [];
    
    if isempty(ind)
       break; 
    end
end

% over-approximate all selected generators that did not get splitted
if ~isempty(ind)

    expMat1 = [expMat1;zeros(1,size(expMat1,2))];
    
    for i = 1:length(ind)
        expMat1(end,ind(i)) = expMat1(index,ind(i));
        expMat1(index,ind(i)) = 0;
    end
end

expMat2 = expMat1;

% add every generator with all-zero exponent matrix to the zonotope center
temp = sum(expMat1,1);
ind = find(temp == 0);
c1 = c1 + sum(G1(:,ind),2);
G1(:,ind) = [];
expMat1(:,ind) = [];

temp = sum(expMat2,1);
ind = find(temp == 0);
c2 = c2 + sum(G2(:,ind),2);
G2(:,ind) = [];
expMat2(:,ind) = [];

% construct the resulting polynomial zonotopes
pZsplit{1} = polyZonotope(c1,G1,pZ.Grest,expMat1);
pZsplit{2} = polyZonotope(c2,G2,pZ.Grest,expMat2);

%------------- END OF CODE --------------