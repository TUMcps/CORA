function [E] = generateRandom(varargin)
% generateRandom - Generates a random ellipsoid
%
% Syntax:  
%    E = generateRandom() Generates an ellipsoid with random orientation,
%        center and dimension
%    E = generateRandom(isdegenerate) Generates a random
%        degenerate/non-degenerate ellipsoid based on given boolean
%        isdegenerate
%    E = generateRandom(isdegenerate,dim) Generates a random
%        degenerate/non-degenerate ellipsoid with dimension dim
%    E = generateRandom(dim,isdegenerate) see above
%
% Inputs:
%   (opt.) isdegenerate - type of ellipsoid
%
% Outputs:
%    E - random ellipsoid
%
% Example: 
%    E = ellipsoid.generateRandom(true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann, Matthias Althoff
% Written:      13-March-2019
% Last update:  02-Sep-2019 (rename generate -> generateRandom)
% Last revision:---
%------------- BEGIN CODE --------------

pd = makedist('Normal','mu',0,'sigma',3); 
isdegenerate=[];
n = ceil(abs(random(pd,1,1))+1.01);

if nargin==1
    if isa(varargin{1},'logical')
        isdegenerate = varargin{1};
    elseif floor(varargin{1})==varargin{1} && varargin{1}>0
        n = varargin{1};
    end
elseif nargin==2
    if isa(varargin{1},'logical') && floor(varargin{2})==varargin{2} && varargin{2}>0
        isdegenerate = varargin{1};
        n = varargin{2};
    elseif isa(varargin{2},'logical') && floor(varargin{1})==varargin{1} && varargin{1}>0
        isdegenerate = varargin{2};
        n = varargin{1};
    else
        error('Wrong input arguments');
    end
elseif nargin>2
    error('Too many input arguments');
end

tmp = random(pd,n,n);
Q = tmp'*tmp;
Q = 1/2*(Q+Q');
if ~isempty(isdegenerate)
    isdegenerate=varargin{1};
    [V,D] = eig(Q);
    [~,ind] = sort(diag(D),'descend');
    D = D(ind,ind);
    V = V(:,ind);
    d = diag(D);
    if isdegenerate && rank(D)==length(D)
        i = ceil((n-1)*rand);
        ind = unique(ceil(n*rand(i,1)));
        d(ind)=0;
        Q = V*diag(d)*V';
    elseif ~isdegenerate && rank(D)<length(D)
        while(1)
            i = length(D)-rank(D);
            d(end-i+1:end) = abs(random(pd,i,1));
            if rank(diag(d))==length(D)
                break;
            end
        end
        Q = V*diag(d)*V';
    end
end

q = random(pd,n,1);
E = ellipsoid(Q,q);
%------------- END OF CODE --------------