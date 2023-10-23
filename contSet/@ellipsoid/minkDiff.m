function E = minkDiff(E,S,varargin)
% minkDiff - computes the Minkowski difference of an ellipsoid as a minuend
%    and a set as a subtrahend
%
% Syntax:
%    E = minkDiff(E,S)
%    E = minkDiff(E,S,mode)
%    E = minkDiff(E,S,L)
%    E = minkDiff(E,S,L,mode)
%
% Inputs:
%    E - ellipsoid object
%    S - set representation/double matrix
%    L - (optional) directions to use for approximation
%    mode - (optional) type of approximation ('inner', 'outer')
%
% Outputs:
%    E - ellipsoid object after Minkowski difference
%
% Example: 
%    E1 = ellipsoid(10*eye(2));
%    E2 = ellipsoid(8*eye(2));
% 
%    res = minkDiff(E1,E2);
% 
%    figure; hold on;
%    plot(E1); plot(E2);
%    plot(res,[1,2],'r');
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

% Authors:       Victor Gassmann
% Written:       09-March-2021
% Last update:   04-July-2022 (VG, class array instead of cell array)
%                09-November-2022 (MW, rename 'minkDiff')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% parsing and checking
% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {S,'att',{'contSet','numeric'}}});

% check dims
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


% subtrahend is the empty set
if all(representsa_(S,'emptySet',eps))
    return;
end

% check arguments
if ~any(strcmp(mode,{'inner','outer'}))
    throw(CORAerror('CORA:wrongValue','fourth',"'inner' or 'outer'"));
end

if size(L,1)~=dim(E)
    throw(CORAerror('CORA:dimensionMismatch',E,L));
end

%% different Minkowski differences
N = length(S);

if isa(S,'double')
    s = sum(S,2);
    E = ellipsoid(E.Q,E.q-s);
    return;
end

if isa(S,'ellipsoid')
    for i=1:N
        E = minkDiffEllipsoid(E,S(i),L,mode);
    end
    return;
end

% ------------------------------ END OF CODE ------------------------------
