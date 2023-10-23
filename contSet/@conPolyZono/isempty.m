function res = isempty(cPZ,varargin)
% isempty - checks if a constrained polynomial zonotope is the empty set
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
%    res - true/false
%
% Example:
%    c = [0;0];
%    G = [1 0 1;0 1 1];
%    E = [1 0 2;0 1 1];
%    A = [1 -1 0; 0 -1 1];
%    b1 = [0; 1]; b2 = [0; 0];
%    EC = [2 0 1; 0 1 0];
%    cPZ1 = conPolyZono(c,G,E,A,b1,EC);
%    cPZ2 = conPolyZono(c,G,E,A,b2,EC);
%
%    res1 = isempty(cPZ1,'linearize',3,7)
%    res2 = isempty(cPZ2,'linearize',3,7)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/isempty, contract

% Authors:       Niklas Kochdumper
% Written:       04-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if independent generators are empty 
if ~isempty(cPZ.GI)
   res = false; return; 
end

% check if constraints exist
if isempty(cPZ.A)
   res = false; return; 
end

% parse input arguments
[method,splits,iter] = setDefaultValues({'linearize',0,1},varargin);

% check input arguments
inputArgsCheck({{cPZ,'att','conPolyZono'};
                {method,'str',{'forwardBackward','linearize',...
                    'polynomial','interval','all'}};
                {splits,'att','numeric',{'scalar','integer'}};
                {iter,'att','numeric',{'scalar','integer'}}});
    
% try to contract the domain to the empty set -> set is empty
temp = ones(length(cPZ.id),1);
dom = interval(-temp,temp);

D = contractPoly(-cPZ.b,cPZ.A,[],cPZ.EC,dom,method,iter,splits);

res = isempty(D);

% ------------------------------ END OF CODE ------------------------------
