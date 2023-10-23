function index = subsetIndex(v,vs)
% subsetIndex - Returns index of each element of vs in v
%
% Syntax:
%    index = subsetIndex(v,vs)
%
% Inputs:
%    v  - vector
%    vs - vector (each component must occur in v)
%
% Outputs:
%    index - index such that v(index) = vs
%
% Example: 
%    v = [5;2;7];
%    vs = [7;5];
%    index = subsetIndex(v,vs) % output is [3;1]
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       25-February-2022 
% Last update:   17-November-2022 (VG, use inputArgsCheck)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({ ...
                {v,'att',{'double'},{'vector'}};
                {vs,'att',{'double'},{'vector'}}...
                });

if ~all(ismember(vs,v))
    throw(CORAerror('CORA:wrongValue','first/second',...
        'At least one component of "v_sub" not in "v"!'));
end

[r,ii_v,~] = intersect(v,vs,'stable');
% since vs is subset of v, r and vs are (up to permutations)
% the same => compute indices to go from r to vs, i.e. find ii_r such that
% r(ii_r) = vs
% the first argument is identical to vs (thanks to 'stable')
[~,~,ii_r] = intersect(vs,r,'stable');
index = ii_v(ii_r);

% ------------------------------ END OF CODE ------------------------------
