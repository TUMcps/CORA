function E = minkDiff(E,S,varargin)
% minkDiff - computes the Minkowski difference of an ellipsoid as a minuend
%    and a set as a subtrahend
%
% Syntax:
%    E = minkDiff(E,S)
%    E = minkDiff(E,S,mode)
%    E = minkDiff(E,S,mode,L)
%
% Inputs:
%    E - ellipsoid object
%    S - set representation/double matrix
%    mode - (optional) type of approximation ('inner', 'outer')
%    L - (optional) directions to use for approximation
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

narginchk(2,4);

[mode,L] = setDefaultValues({'outer',zeros(dim(E),0)},varargin);

%% parsing and checking
% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {S,'att',{'contSet','numeric','cell'}}; ...
                {mode,'str',{'inner','outer'}}; ...
                {L,'att','numeric'}});

% check dims
equalDimCheck(E,S);
equalDimCheck(E,L);

% subtrahend is the empty set
if ~iscell(S) && representsa_(S,'emptySet',eps)
    E = fullspace(dim(E));
    return;
end

% different Minkowski differences
if isnumeric(S) && iscolumn(S)
    s = sum(S,2);
    E = ellipsoid(E.Q,E.q-s);
    return;
end

% ensure that all subtrahends are ellipsoids
if iscell(S) && all(cellfun(@(S_i) isa(S_i,'ellipsoid'),S,'UniformOutput',true))
    for i=1:length(S)
        E = priv_minkDiffEllipsoid(E,S{i},L,mode);
    end
    return
end

if isa(S,'ellipsoid')
    E = priv_minkDiffEllipsoid(E,S,L,mode);
    return;
end

throw(CORAerror('CORA:noops',E,S));

% ------------------------------ END OF CODE ------------------------------
