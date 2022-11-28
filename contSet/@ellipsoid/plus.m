function E = plus(E,S,varargin)
% plus - Overloaded '+' operator for approximating the Minkowski sum of an
%    ellipsoid and another set 
%
% Syntax:  
%    E = plus(E,S)
%    E = plus(E,S,mode)
%    E = plus(E,S,L)
%    E = plus(E,S,L,mode)
%
% Inputs:
%    E - ellipsoid object
%    S - set representation/double matrix
%    L - (optional) directions to use for approximation
%    mode - (optional) type of approximation
%               'inner'
%               'outer':
%               'outer:halder': available when L is empty
%
% Outputs:
%    E - ellipsoid object after Minkowski sum
%
% Example: 
%    E1 = ellipsoid(eye(2),[1;-1]);
%    E2 = ellipsoid(diag([1,2]));
%    Ep = E1 + E2;
%    figure; hold on
%    plot(E1); plot(E2);
%    plot(Ep,[1,2],'r');
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal
%       toolbox (ET). In Proceedings of the 45th IEEE Conference on
%       Decision and Control (pp. 1498-1503). IEEE.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      09-March-2021
% Last update:  04-July-2022 (VG: class array instead of cell array)
% Last revision:---

%------------- BEGIN CODE --------------

%% parsing & checking
% make sure first argument is class argument 
[E,S] = findClassArg(E,S,'ellipsoid');

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {S,'att',{'contSet','numeric'}}});

% Minkowski addition with empty set
if isempty(E)
    return;
elseif (isnumeric(S) && isempty(S)) || (isa(S,'contSet') && any(isemptyobject(S)))
    E = ellipsoid(); return
end

% dimension check
equalDimCheck(E,S);

% parse arguments
if isempty(varargin)
    L = zeros(dim(E),0);
    mode = 'outer';
elseif length(varargin)==1
    if isa(varargin{1},'char')
        mode = varargin{1};
        L = zeros(dim(E),0);
    elseif isa(varargin{1},'double')
        L = varargin{1};
        mode = 'outer';
    else
        throw(CORAerror('CORA:wrongValue','third',"be of type 'double' or 'char'"));
    end
elseif length(varargin)==2
    if ~isa(varargin{1},'double')
        throw(CORAerror('CORA:wrongValue','third',"be of type 'double'"));
    end
    if ~isa(varargin{2},'char')
        throw(CORAerror('CORA:wrongValue','fourth',"be of type 'char'"));
    end
    L = varargin{1};
    mode = varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',4));
end


if all(isempty(S))
    return;
end

% check arguments
if ~strcmp(mode,'outer') && ~strcmp(mode,'inner')
    if length(varargin)==2
        throw(CORAerror('CORA:wrongValue','fourth',"'inner' or 'outer'"));
    else
        if ~isa(S,'ellipsoid') || ~strcmp(mode,'outer:halder')
            throw(CORAerror('CORA:wrongValue','third',"'inner' or 'outer'"));
        end
    end
end

if size(L,1)~=dim(E)
    throw(CORAerror('CORA:dimensionMismatch',E,L));
end

N = length(S);

%% different Minkowski additions
if isa(S,'double')
    s = sum(S,2);
    E = ellipsoid(E.Q,E.q+s);
    return;
end

if isa(S,'conPolyZono')
    E = S(1) + E; 
    for i=2:N
        E = S(i) + E; 
    end
    return; 
end

if isa(S,'ellipsoid')
   E = plusEllipsoid([E;S(:)],L,mode);
   return;
end

% throw error for all other combinations
throw(CORAerror('CORA:noops',E,S));

%------------- END OF CODE --------------