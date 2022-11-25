function I = interval(pZ,varargin)
% interval - Over-approximates a polynomial Zonotope by an interval
%
% Syntax:  
%    I = interval(pZ)
%    I = interval(pZ,method)
%
% Inputs:
%    qZ - polyZonotope object
%    method - method used to calculate the bounds for all dimensions
%              'interval': interval arithmetic
%              'split': split set multiple times
%              'bnb': taylor models with "branch and bound" algorithm
%              'bnbAdv': taylor models with advandced bnb-algorithm
%              'globOpt': verified global optimization 
%              'bernstein': conversion to a bernstein polynomial
%
% Outputs:
%    I - interval object
%
% Example: 
%    pZ = polyZonotope([0;0],[2 -1 2;0 -2 -3],[],[1 0 1;0 1 3]);
%
%    int1 = interval(pZ);
%    int2 = interval(pZ,'bernstein');
%
%    figure
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(int1,[1,2],'b');
%    xlim([-6,6]);
%    ylim([-6,6]);
%
%    figure
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(int2,[1,2],'g');
%    xlim([-6,6]);
%    ylim([-6,6]);
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, supportFunc

% Author:       Niklas Kochdumper
% Written:      23-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------p
% parse input arguemnts
method = 'interval';
if nargin == 2
   method = varargin{1}; 
end

% compute over-approximating interval with the selected method
if strcmp(method,'interval')
    
    I = interval(zonotope(pZ));
    
elseif strcmp(method,'bernstein')
    
    p = size(pZ.expMat,1);    
    dom = interval(-ones(p,1),ones(p,1));  
    
    % dependent generators: convert to bernstein polynomial
    B = poly2bernstein(pZ.G,pZ.expMat,dom);
       
    infi = cellfun(@(x) min(min(x)),B);
    sup = cellfun(@(x) max(max(x)),B);
    
    int1 = interval(infi,sup);
    
    % independent generators: enclose zonotope with interval
    int2 = interval(zonotope([pZ.c,pZ.Grest]));
    
    I = int1 + int2;
    
else
   
   % initialize variables 
   n = length(pZ.c);
   e = zeros(n,1);
   I = interval(e,e);
   
   % loop over all system dimensions
   for i = 1:n
      e_ = e;
      e_(i) = 1;
      I(i) = supportFunc(pZ,e_,'range',method);
   end
end

%------------- END OF CODE --------------