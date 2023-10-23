function Z = zonotope(cPZ,varargin)
% zonotope - Computes a zonotope that over-approximates the constrained
%    polynomial zonotope
%
% Syntax:
%    Z = zonotope(cPZ)
%    Z = zonotope(cPZ,method)
%
% Inputs:
%    cPZ - conPolyZono object
%    method - algorithm used for contraction ('forwardBackward',
%             'linearize', 'polynomial', 'interval', 'all', or 'none')
%
% Outputs:
%    Z - zonotope object
%
% Example:  
%    c = [0;0];
%    G = [2 1 2 1; 0 2 2 1];
%    E = [1 0 2 0; 0 1 1 0; 0 0 0 1];
%    A = [1 -0.5 0.5];
%    b = 0.5;
%    EC = [1 0 0; 0 1 2; 0 1 0];
%
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    Z = zonotope(cPZ);
%   
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','r','Splits',10);
%    plot(Z);
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope, interval, polyZonotope

% Authors:       Niklas Kochdumper
% Written:       07-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
method = setDefaultValues({'linearize'},varargin);

% check input arguments
inputArgsCheck({{cPZ,'att','conPolyZono'};
                {method,'str',{'forwardBackward','linearize',...
                    'polynomial','interval','all','none'}}});

% contract the domain for the factors based on polynomial constraints
temp = ones(length(cPZ.id),1);
dom = interval(-temp,temp);

if ~isempty(cPZ.A) && ~strcmp(method,'none')
	dom = contractPoly(-cPZ.b,cPZ.A,[],cPZ.EC,dom,method);
end

if representsa_(dom,'emptySet',eps)
    throw(CORAerror('CORA:emptySet'));
end

% construct enclosing zonotope
pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.GI,cPZ.E,cPZ.id);
S = getSubset(pZ,pZ.id,dom);
Z = zonotope(S);
    
% ------------------------------ END OF CODE ------------------------------
