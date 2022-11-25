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
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%
%    figure
%    hold on
%    plot(pZsplit{1},[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(pZsplit{2},[1,2],'g','Filled',true,'EdgeColor','none');
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: split, splitLongestGen

% Author:       Niklas Kochdumper
% Written:      24-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% find selected dependent factor
index = find(pZ.id == ind);
if isempty(index)
   error('Dependent factor with the specified identifier does not exist!'); 
end

% determine all generators in which the selected dependent factor occurs
ind = find(pZ.expMat(index,:) > 0);

% parse input arguments
if nargin == 3
   polyOrd = varargin{1}; 
else
   polyOrd = max(pZ.expMat(index,ind));
end

% determine coefficients of polynomials that correspond to (0.5 + 0.5x)^p
polyCoeff1 = cell(max(2,polyOrd),1);
polyCoeff1{1} = 0.5 * [1, 1];
polyCoeff1{2} = 0.25 * [1, 2, 1];

if polyOrd >= 3
   syms x
   y = (1+x);
end

for i = 3:polyOrd
    temp = y^i;
    polyCoeff1{i} = 0.5^i * eval(coeffs(temp));
end


% determine coefficients of polynomials that correspond to (-0.5 + 0.5x)^p
polyCoeff2 = cell(max(2,polyOrd),1);
polyCoeff2{1} = 0.5 * [-1, 1];
polyCoeff2{2} = 0.25 * [1, -2, 1];

if polyOrd >= 3
   syms x
   y = (-1+x);
end

for i = 3:polyOrd
    temp = y^i;
    polyCoeff2{i} = 0.5^i * eval(coeffs(temp));
end

% construct the modified generators for the splitted zonotopes
expMat1 = pZ.expMat;
expMat2 = pZ.expMat;
G1 = pZ.G;
G2 = pZ.G;
c1 = pZ.c;
c2 = pZ.c;

for i = 1:polyOrd
   
    indTemp = find(pZ.expMat(index,ind) == i);
    
    for j = 1:length(indTemp)
       
        g = pZ.G(:,ind(indTemp(j)));
        e = pZ.expMat(:,ind(indTemp(j)));
       
        % first splitted zonotope
        coef = polyCoeff1{i};
        
        G1(:,ind(indTemp(j))) = coef(1) * g;
        expMat1(index,ind(indTemp(j))) = 0;
        
        G1 = [G1, g*coef(2:end)];
        E_ = e* ones(1,length(coef)-1);
        E_(index,:) = 1:length(coef)-1;
        expMat1 = [expMat1, E_];
        
        
        % second splitted zonotope
        coef = polyCoeff2{i};
        
        G2(:,ind(indTemp(j))) = coef(1) * g;
        expMat2(index,ind(indTemp(j))) = 0;
        
        G2 = [G2, g*coef(2:end)];
        E_ = e* ones(1,length(coef)-1);
        E_(index,:) = 1:length(coef)-1;
        expMat2 = [expMat2, E_];
        
    end
    
    % remove the finished indizes from the list
    ind(indTemp) = [];
    
    if isempty(ind)
       break; 
    end
end

% over-approximate all selected generators that did not get splitted
if ~isempty(ind)

    expMat1 = [expMat1;zeros(1,size(expMat1,2))];
    expMat2 = [expMat2;zeros(1,size(expMat2,2))];
    
    for i = 1:length(ind)

        expMat1(end,ind(i)) = expMat1(index,ind(i));
        expMat1(index,ind(i)) = 0;
        
        expMat2(end,ind(i)) = expMat2(index,ind(i));
        expMat2(index,ind(i)) = 0;  
    end
end

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

% remove redundant exponents
[expMat1,G1] = removeRedundantExponents(expMat1,G1);
[expMat2,G2] = removeRedundantExponents(expMat2,G2);

% construct the resulting polynomial zonotopes
pZsplit{1} = polyZonotope(c1,G1,pZ.Grest,expMat1);
pZsplit{2} = polyZonotope(c2,G2,pZ.Grest,expMat2);


%------------- END OF CODE --------------