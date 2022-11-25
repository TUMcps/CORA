function res = isempty(cPZ,varargin)
% isempty - check if a constrained polynomial zonotope is empty
%
% Syntax:  
%    res = isempty(cPZ)
%    res = isempty(cPZ,method)
%    res = isempty(cPZ,method,iter)
%    res = isempty(cPZ,method,iter,splits)
%
% Inputs:
%    cPZ - conPolyZono object
%    method - algorithm used for contraction ('forwardBackward',
%            'linearize', 'polynomial', 'interval', or 'all')
%    iter - number of iteration (integer > 0 or 'fixpoint')
%    splits - number of recursive splits (integer > 0)
%
% Outputs:
%   res - true is set is empty, false otherwise
%
% Example:
%    c = [0;0];
%    G = [1 0 1;0 1 1];
%    expMat = [1 0 2;0 1 1];
%    A = [1 -1 0; 0 -1 1];
%    b1 = [0; 1]; b2 = [0; 0];
%    expMat_ = [2 0 1; 0 1 0];
%    cPZ1 = conPolyZono(c,G,expMat,A,b1,expMat_);
%    cPZ2 = conPolyZono(c,G,expMat,A,b2,expMat_);
%
%    res1 = isempty(cPZ1,'linearize',3,7)
%    res2 = isempty(cPZ2,'linearize',3,7)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isempty, contract

% Author:       Niklas Kochdumper
% Written:      04-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % check if independent generators are empty 
    if ~isempty(cPZ.Grest)
       res = false; return; 
    end
    
    % check if constraints exist
    if isempty(cPZ.A)
       res = false; return; 
    end
    
    % parse input arguments
    method = 'linearize';
    splits = 0; 
    iter = 1;
    
    if nargin > 1 && ~isempty(varargin{1})
       method = varargin{1}; 
    end
    if nargin > 2 && ~isempty(varargin{2})
       iter = varargin{2};
    end
    if nargin > 3 && ~isempty(varargin{3})
       splits = varargin{3};
    end
    
    % try to contract the domain to the empty set -> set is empty
    temp = ones(length(cPZ.id),1);
    dom = interval(-temp,temp);
    
    D = contractPoly(-cPZ.b,cPZ.A,[],cPZ.expMat_,dom,method,iter,splits);
    
    res = isempty(D);
end

%------------- END OF CODE --------------