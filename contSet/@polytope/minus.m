function P_out = minus(P,S,varargin)
% minus - dummy function to alert users of the difference in meaning
%    between 'minus' for range bounding and 'minkDiff' for the Minkowski
%    difference; for numerical vectors as subtrahends, these operations are
%    equivalent, so we compute it here nonetheless
%
% Syntax:
%    P_out = minus(P,S)
%
% Inputs:
%    P - polytope object
%    S - numeric
%
% Outputs:
%    P_out - polytope object
%
% Example: 
%    P = polytope([1 0; -1 1; -1 -1],[1;1;1]);
%    v = [1;0];
%    P_ = P - v;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minkDiff

% Authors:       Mark Wetzlinger
% Written:       09-November-2022
% Last update:   03-October-2024 (MW, use private method)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only supported for subtrahends that are numeric vectors
if isnumeric(S) && iscolumn(S)
    % skip dimension check for speed... copy polytope:
    P_out = polytope(P);
    if P.isHRep.val
        [A,b,Ae,be] = priv_plus_minus_vector(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val,-S);
        P_out.A_.val = A; P_out.b_.val = b;
        P_out.Ae_.val = Ae; P_out.be_.val = be;
    end
    if P.isVRep.val
        V = P.V_.val - S;
        P_out.V_.val = V;
    end

    % copy properties
    P_out = priv_copyProperties(P,P_out,'noV');

else
    % throw error
    throw(CORAerror('CORA:notSupported',...
        ['The function ''minus'' is not implemented for the class polytope except for vectors as a subtrahend.\n', ...
        'If you require to compute the Minkowski difference, use ''minkDiff'' instead.']));
end

% ------------------------------ END OF CODE ------------------------------
