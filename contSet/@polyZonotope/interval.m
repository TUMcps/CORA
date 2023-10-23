function I = interval(pZ,varargin)
% interval - Over-approximates a polynomial zonotope by an interval
%
% Syntax:
%    I = interval(pZ)
%    I = interval(pZ,method)
%
% Inputs:
%    pZ - polyZonotope object
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
%    I1 = interval(pZ);
%    I2 = interval(pZ,'bernstein');
%
%    figure; hold on; xlim([-6,6]); ylim([-6,6]);
%    plot(pZ,[1,2],'FaceColor','r');
%    plot(I1,[1,2],'b');
%
%    figure; hold on; xlim([-6,6]); ylim([-6,6]);
%    plot(pZ,[1,2],'FaceColor','r');
%    plot(I2,[1,2],'g');
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, supportFunc

% Authors:       Niklas Kochdumper
% Written:       23-March-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguemnts
method = setDefaultValues({'interval'},varargin);

% check input arguments
inputArgsCheck({{pZ,'att','polyZonotope'};
                {method,'str',{'interval','split','bnb','bnbAdv',...
                    'globOpt','bernstein'}}});

% compute over-approximating interval with the selected method
if strcmp(method,'interval')
    
    I = interval(zonotope(pZ));
    
elseif strcmp(method,'bernstein')
    
    p = size(pZ.E,1);    
    dom = interval(-ones(p,1),ones(p,1));  
    
    % dependent generators: convert to bernstein polynomial
    B = poly2bernstein(pZ.G,pZ.E,dom);
       
    infi = cellfun(@(x) min(min(x)),B);
    sup = cellfun(@(x) max(max(x)),B);
    
    I1 = interval(infi,sup);
    
    % independent generators: enclose zonotope with interval
    I2 = interval(zonotope([pZ.c,pZ.GI]));
    
    I = I1 + I2;
    
else
   
    % initialize variables 
    n = length(pZ.c);
    e = zeros(n,1);
    I = interval(e,e);
    
    % loop over all system dimensions
    for i = 1:n
        e_ = e;
        e_(i) = 1;
        I(i) = supportFunc_(pZ,e_,'range',method,8,1e-3);
    end

end

% ------------------------------ END OF CODE ------------------------------
