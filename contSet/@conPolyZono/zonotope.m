function res = zonotope(obj,varargin)
% zonotope - Computes a zonotope that over-approximates the conPolyZonotope
%            object
%
% Syntax:  
%    res = zonotope(obj)
%    res = zonotope(obj,method)
%
% Inputs:
%    obj - conPolyZono object
%    method - algorithm used for contraction ('forwardBackward',
%             'linearize', 'polynomial', 'interval', 'all', or 'none')
%
% Outputs:
%    res - zonotope object
%
% Example:  
%    c = [0;0];
%    G = [2 1 2 1; 0 2 2 1];
%    expMat = [1 0 2 0; 0 1 1 0; 0 0 0 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    expMat_ = [1 0 0; 0 1 2; 0 1 0];
%
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    Z = zonotope(cPZ);
%   
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Splits',10,'Filled',true,'EdgeColor','none');
%    plot(Z);
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope, interval, polyZonotope

% Author:       Niklas Kochdumper
% Written:      07-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    method = 'linearize';
    
    if nargin > 1 && ~isempty(varargin{1})
       method = varargin{1}; 
    end

    % contract the domain for the factors based on polynomial constraints
    temp = ones(length(obj.id),1);
    dom = interval(-temp,temp);
    
    if ~isempty(obj.A) && ~strcmp(method,'none')
    	dom = contractPoly(-obj.b,obj.A,[],obj.expMat_,dom,method);
    end
    
    if isempty(dom)
       [msg,id] = errEmptySet();
       error(id,msg); 
    end
    
    % construct enclosing zonotope
    pZ = polyZonotope(obj.c,obj.G,obj.Grest,obj.expMat,obj.id);
    
    S = getSubset(pZ,pZ.id,dom);
    
    res = zonotope(S);
end
    
%------------- END OF CODE --------------