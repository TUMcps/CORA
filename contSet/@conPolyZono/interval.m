function res = interval(obj,varargin)
% interval - computes an enclosing interval of a constrained polynomial
%            zonotope
%
% Syntax:  
%    res = interval(obj)
%    res = interval(obj,method)
%
% Inputs:
%    obj - conPolyZono object
%    method - method that is used to calculate the interval enclosure
%              'conZonotope': conversion to a constrained zonotope
%              'interval': interval arithmetic
%              'split': split set multiple times
%              'quadProg': quadratic programming
%
% Outputs:
%    res - interval object
%
% Example: 
%    A = 1/8 * [-10 2 2 3 3];
%    b = -3/8;
%    expMat_ = [1 0 1 2 0; 0 1 1 0 2; 0 0 0 0 0];
%    c = [0;0];
%    G = [1 0 1 -1/4;0 1 -1 1/4];
%    expMat = [1 0 2 0;0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    int = interval(cPZ);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Filled',true,'EdgeColor','none','Splits',20);
%    plot(int,[1,2],'b');
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: supportFunc, conZonotope

% Author:       Niklas Kochdumper
% Written:      14-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    method = 'conZonotope';
    
    if nargin > 1
       method = varargin{1}; 
    end
    
    % compute enclosing interval with the spe
    if strcmp(method,'conZonotope')
        
        % compute conZonotope enclosure of conPolyZono object
        zono = conZonotope(obj);
        
        % enclose zonotope with an interval
        res = interval(zono);
        
    elseif strcmp(method,'interval')
        
        % compute zonotope enclosure of conPolyZono object
        zono = zonotope(obj);
        
        % enclose zonotope with an interval
        res = interval(zono);
        
    elseif strcmp(method,'split') || strcmp(method,'quadProg')
        
        n = dim(obj);
        res = interval(zeros(n,1));

        % loop over all dimensions
        for i = 1:n
            
            % construct unit vector
            temp = zeros(n,1);
            temp(i) = 1;

            % calculate bounds
            res(i) = supportFunc(obj,temp,'range',method);
        end
        
    else
        error('Wrong value for input argument "method"!'); 
    end
end

%------------- END OF CODE --------------