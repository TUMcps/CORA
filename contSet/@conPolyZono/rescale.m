function res = rescale(obj,varargin)
% rescale - rescale a constrained polynomial zonotope object by computing 
%           a shrinked factor domain resulting from the constraints
%
% Syntax:  
%    res = rescale(obj)
%    res = rescale(obj,method)
%
% Inputs:
%    obj - conPolyZono object
%    method - algorithm used for contraction ('forwardBackward',
%             'linearize', 'polynomial', 'interval', or 'all')
%
% Example:
%    c = [0;0];
%    G = [1 0 0 0.2;0 -2 1 0.2];
%    expMat = [1 1 2 0;0 1 1 0;0 0 0 1];
%    A = [1 2 1 2 0.75];
%    b = -1.75;
%    expMat_ = [2 1 0 0 0;0 0 2 1 0;0 0 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    cPZ_ = rescale(cPZ,'forwardBackward');
%
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Filled',true,'EdgeColor','none','Splits',20);
%    plot(polyZonotope(c,G,[],expMat),[1,2],'b','Splits',15);
%    plot(polyZonotope(cPZ_.c,cPZ_.G,[],cPZ_.expMat),[1,2],'g','Splits',15);
%
% Outputs:
%    res - rescaled conPolyZono object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/rescale

% Author:       Niklas Kochdumper
% Written:      06-November-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE -------------

    % parse input arguments
    method = 'forwardBackward';
    
    if nargin > 1
        method = varargin{1};
    end
    
    % contract the domain for the factors based on polynomial constraints
    temp = ones(length(obj.id),1);
    dom = interval(-temp,temp);
    
    if ~isempty(obj.A)
    	dom = contractPoly(-obj.b,obj.A,[],obj.expMat_,dom,method);
    end
    
    if isempty(dom)
       [msg,id] = errEmptySet();
       error(id,msg); 
    end
    
    % reexpand the conPolyZono object so that the range for the factors is
    % again between [-1,1]
    res = getSubset(obj,obj.id,dom);
    
end
    
%------------- END OF CODE --------------