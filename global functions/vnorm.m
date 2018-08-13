function y = vnorm(A,varargin)
% VNORM - Return the vector norm along specified dimension of A
%
%   VNORM(A) returns the 2-norm along the first non-singleton
%   dimension of A
%   VNORM(A,dim) return the 2-norm along the dimension 'dim'
%   VNORM(A,dim,normtype) returns the norm specified by normtype
%   along the dimension 'dim'
%   VNORM(A,[],normtype) returns the norm specified by normtype along
%   the first non-singleton dimension of A
% 
%   normtype may be one of {inf,-inf,positive integer}.
%   For a given vector, v, these norms are defined as
%   inf: max(abs(v))
%   -inf: min(abs(v))
%   p (where p is a positive integer): sum(abs(v).^p)^(1/p)
%   
%   Examples:
%       A =  [8 1 6; 3 5 7; 4 -9 2];
%
%       %Columnwise 2-norm (Euclidean norm)
%       vnorm(A,1) = [9.4340 10.3441 9.4340];
%       vnorm(A,[],2) % Same as above (since first non-singleton dimensions
%                     % is columnwise and default norm is 2-norm.
%       vnorm(A,[],[])% Again, same as above
%
%       % Row-wise maximum of absolute values
%       vnorm(A,2,inf) = [8 7 9]'; 
%
%       % Columnwise minimum of absolute values
%       vnorm(A,[],-inf) = [3 1 2];
%       
%       % Error: Use the inf type and not the string 'inf'
%       vnorm(A,[],'inf')   % Wrong
%       vnorm(A,[],inf)     % Correct

dim = [];
ntype = [];

if nargin>1
    dim = varargin{1};
    if isempty(dim)
       idx = find(size(A)~=1);
       dim = idx(1);
    elseif dim~=floor(dim) || dim<1
        error('Dimension must be positive integer');
    end
    if nargin>2
        ntype = varargin{2};
    end
end


if isempty(ntype)
    y = sqrt(sum( abs(A).^2 , dim) );
elseif ntype==1
    y = sum( abs(A) , dim );    
elseif isinf(ntype)
    if ntype > 0
        y=max(abs(A), [], dim);
    else
        y=min(abs(A), [], dim);
    end
elseif ntype~=floor(ntype) || ntype<1
    error(['Norm type must be one of inf,-inf or a positive ' ...
           'integer']);
else 
    y = (sum( abs(A).^ntype , dim) ).^(1/ntype);
end 
