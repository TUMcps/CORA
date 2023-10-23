function res = rescale(cPZ,varargin)
% rescale - rescale a constrained polynomial zonotope object by computing 
%    a shrinked factor domain resulting from the constraints
%
% Syntax:
%    res = rescale(cPZ)
%    res = rescale(cPZ,method)
%
% Inputs:
%    cPZ - conPolyZono object
%    method - algorithm used for contraction ('forwardBackward',
%             'linearize', 'polynomial', 'interval', or 'all')
%
% Example:
%    c = [0;0];
%    G = [1 0 0 0.2;0 -2 1 0.2];
%    E = [1 1 2 0;0 1 1 0;0 0 0 1];
%    A = [1 2 1 2 0.75];
%    b = -1.75;
%    EC = [2 1 0 0 0;0 0 2 1 0;0 0 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    cPZ_ = rescale(cPZ,'forwardBackward');
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%    plot(polyZonotope(c,G,[],E),[1,2],'b','Splits',12);
%    plot(polyZonotope(cPZ_.c,cPZ_.G,[],cPZ_.E),[1,2],'g','Splits',12);
%
% Outputs:
%    res - rescaled conPolyZono object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/rescale

% Authors:       Niklas Kochdumper
% Written:       06-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
method = setDefaultValues({'forwardBackward'},varargin);

% check input arguments
inputArgsCheck({{cPZ,'att','conPolyZono'};
                {method,'str',{'forwardBackward','linearize',...
                    'polynomial','interval','all'}}});    

% contract the domain for the factors based on polynomial constraints
temp = ones(length(cPZ.id),1);
dom = interval(-temp,temp);

if ~isempty(cPZ.A)
	dom = contractPoly(-cPZ.b,cPZ.A,[],cPZ.EC,dom,method);
end

if representsa_(dom,'emptySet',eps)
    throw(CORAerror('CORA:emptySet'));
end

% reexpand the conPolyZono object so that the range for the factors is
% again between [-1,1]
res = getSubset(cPZ,cPZ.id,dom);
    
% ------------------------------ END OF CODE ------------------------------
