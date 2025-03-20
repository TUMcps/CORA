function [val,x] = supportFunc_(P,dir,type,varargin)
% supportFunc_ - Calculate the upper or lower bound of a polytope object
%    along a certain direction
%
% Syntax:
%    val = supportFunc_(P,dir,type)
%    [val,x] = supportFunc_(P,dir,type)
%
% Inputs:
%    P - polytope object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower', 'upper', 'range')
%
% Outputs:
%    val - bound of the constraind zonotope in the specified direction
%    x - support vector
% 
% Example:
%    A = [1 0; -1 1; -1 -1]; b = [1;1;1];
%    P = polytope(A,b);
%    [val,x] = supportFunc(P,[0;-1],'upper');
% 
% Reference:
%    [1] M. Wetzlinger, V. Kotsev, A. Kulmburg, M. Althoff. "Implementation
%        of Polyhedral Operations in CORA 2024", ARCH'24.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, conZonotope/supportFunc_

% Authors:       Niklas Kochdumper, Victor Gassmann, Mark Wetzlinger
% Written:       19-November-2019
% Last update:   16-March-2021 (added unbounded support)
%                10-December-2022 (MW, add 'range')
%                13-December-2022 (MW, add call to MOSEK)
%                15-November-2023 (MW, computation for vertex representation)
%                09-July-2024 (TL, added fallback during linprog)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check fullspace
if representsa_(P,'fullspace',0)
    [val,x] = aux_supportFunc_fullspace(dir,type);
    return
end

% check if vertex representation given (skips linear program)
if P.isVRep.val
    [val,x] = aux_supportFunc_V(P,dir,type);
    return
end

% compute support function and support vector using linear program
if strcmp(type,'lower') || strcmp(type,'upper')
    [val,x] = aux_solveLinProg(P,dir,type);
elseif strcmp(type,'range')
    [val_upper,x_upper] = aux_solveLinProg(P,dir,'upper');
    [val_lower,x_lower] = aux_solveLinProg(P,dir,'lower');
    % combine values
    val = interval(val_lower,val_upper);
    x = [x_lower x_upper];
end

end


% Auxiliary functions -----------------------------------------------------

function [val,x] = aux_supportFunc_fullspace(dir,type)
% support function of a polytope that represents R^n

x = Inf .* sign(dir);
x(isnan(x)) = 0;
switch type
    case 'upper'
        val = Inf;
    case 'lower'
        val = -Inf;
    case 'range'
        val = interval(-Inf,Inf);
end

end

function [val,x] = aux_supportFunc_V(P,dir,type)
% computation of the support function and support vector of a V-polytope
% according to [1, (32)] and [1, (34)], respectively

if strcmp(type,'upper')
    vals = dir' * P.V_.val;
    [val,idx] = max(vals);
    x = P.V_.val(:,idx);
elseif strcmp(type,'lower')
    vals = -dir' * P.V_.val;
    [val,idx] = min(-vals);
    x = P.V_.val(:,idx);
elseif strcmp(type,'range')
    % compute upper and lower bound
    vals_upper = dir' * P.V_.val;
    vals_lower = -dir' * P.V_.val;
    [val_upper,idx_upper] = max(vals_upper);
    [val_lower,idx_lower] = min(-vals_lower);
    val = interval(val_lower, val_upper);
    x = [P.V_.val(:,idx_lower), P.V_.val(:,idx_upper)];
end

end

function [val,x] = aux_solveLinProg(P,dir,type)
% support function evaluated according to [1, (33)]
% support vector computed according to [1, (35)]

try
    [val,x] = priv_supportFunc(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val,dir,type);
catch ME
    if ~startsWith(ME.identifier,'CORA:solverIssue')
        rethrow(ME);
    end

    % normalize direction with magnitude of constraints 
    % for numeric stability
    norm = max(abs([P.A_.val P.b_.val;P.Ae_.val P.be_.val]),[],'all');
    dir = dir / norm;

    % try to solve linear program
    [val,x] = priv_supportFunc(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val,dir,type);

    % correct val
    val = val * sum((dir~=0) * norm);
end

end

% ------------------------------ END OF CODE ------------------------------
