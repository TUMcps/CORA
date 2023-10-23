function I = interval(cPZ,varargin)
% interval - encloses a constrained polynomial zonotope by an interval
%
% Syntax:
%    I = interval(cPZ)
%    I = interval(cPZ,method)
%
% Inputs:
%    cPZ - conPolyZono object
%    method - method that is used to calculate the interval enclosure
%              'conZonotope': conversion to a constrained zonotope
%              'interval': interval arithmetic
%              'split': split set multiple times
%              'quadProg': quadratic programming
%
% Outputs:
%    I - interval object
%
% Example: 
%    A = 1/8 * [-10 2 2 3 3];
%    b = -3/8;
%    EC = [1 0 1 2 0; 0 1 1 0 2; 0 0 0 0 0];
%    c = [0;0];
%    G = [1 0 1 -1/4;0 1 -1 1/4];
%    E = [1 0 2 0;0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    int = interval(cPZ);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Splits',12);
%    plot(int,[1,2],'b');
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: supportFunc, conZonotope

% Authors:       Niklas Kochdumper
% Written:       14-August-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
method = setDefaultValues({'conZonotope'},varargin);  

% check input arguments
inputArgsCheck({{cPZ,'att','conPolyZono'};
                {method,'str',{'conZonotope','interval','split','quadProg'}}});

% compute enclosing interval with the spe
if strcmp(method,'conZonotope')
    
    % compute conZonotope enclosure of conPolyZono object
    zono = conZonotope(cPZ);
    
    % enclose zonotope with an interval
    I = interval(zono);
    
elseif strcmp(method,'interval')
    
    % compute zonotope enclosure of conPolyZono object
    zono = zonotope(cPZ);
    
    % enclose zonotope with an interval
    I = interval(zono);
    
elseif strcmp(method,'split') || strcmp(method,'quadProg')
    
    n = dim(cPZ);
    I = interval(zeros(n,1));

    % loop over all dimensions
    for i = 1:n
        
        % construct unit vector
        temp = zeros(n,1);
        temp(i) = 1;

        % calculate bounds
        I(i) = supportFunc_(cPZ,temp,'range',method,8);
    end
    
end

% ------------------------------ END OF CODE ------------------------------
