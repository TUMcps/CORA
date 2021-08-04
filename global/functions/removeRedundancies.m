function varargout = removeRedundancies(varargin)
if nargin==1
    if ~isa(varargin{1},'polyZonotope')
        error('polyZonotope required');
    end
    pZ = varargin{1};
    expMat = pZ.expMat;
    G = pZ.G;
    % also remove redundancies in Grest
    ind_rest = all(pZ.Grest==0,1);
    Grest = pZ.Grest(:,~ind_rest);
elseif nargin==2
    expMat = varargin{1};
    G = varargin{2};
    if ~isa(expMat,'double') || (~isa(G,'double') && ~isa(G,'sym'))
        error('Wrong argument type(s)');
    end
    if size(expMat,2)~=size(G,2)
        error('G and expMat  have to have same number of columns');
    end
else
    error('Wrong number of arguments');
end
if ~isempty(G)
    n = size(G,1);
    pZtmp = polyZonotope(zeros(n,1),G,zeros(n,1),expMat);
    % to access removeRedundantExponents
    expMat = pZtmp.expMat;
    G = pZtmp.G;
    ind = all(G==0,1);
    G(:,ind) = [];
    expMat(:,ind) = [];
else
    G = zeros(length(pZ.c),0);
    expMat = zeros(size(expMat,1),0);
end
if nargin==1
    varargout{1} = polyZonotope(pZ.c,G,Grest,expMat,pZ.id);
else
    varargout{1} = expMat;
    varargout{2} = G;
end
