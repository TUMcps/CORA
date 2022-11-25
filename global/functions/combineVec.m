function Y = combineVec(varargin)
% returns all possible (cartesian product) combinations of arguments
% (same behaviour as of combvec)
    
if ~all(cellfun(@(y)isa(y,'double'),varargin))
    error('All input arguments need to be of type "double"!');
end

if nargin==0
    Y = [];
    return;
end
if nargin==1
    Y = varargin{1};
    return;
end

if nargin>2
    Y = combineVec(varargin{1},combineVec(varargin{2},varargin{3}));
    return;
end
Y1 = varargin{1};
Y2 = varargin{2};
[n1,N1] = size(Y1);
[n2,N2] = size(Y2);

Y = zeros(n1+n2,N1*N2);
for i=1:N2
    for j=1:N1
        Y(:,(i-1)*N1 + j) = [Y1(:,j);Y2(:,i)];
    end
end