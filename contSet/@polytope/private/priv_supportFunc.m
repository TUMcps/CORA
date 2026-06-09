function [val,x] = priv_supportFunc(A,b,Ae,be,dir,type)
% priv_supportFunc - computes the support function of a polytope defined by
%    inequality and equality constraints via linear programming
%
% Syntax:
%    [val,x] = priv_supportFunc(A,b,Ae,be,dir,type)
%
% Inputs:
%    A - inequality constraint matrix
%    b - inequality constraint offset
%    Ae - equality constraint matrix
%    be - equality constraint offset
%    dir - direction
%    type - 'upper', 'lower', or 'range'
%
% Outputs:
%    val - value of the support function (interval for 'range')
%    x - support vector (matrix [x_lower x_upper] for 'range')
%
% Other m-files required: supportFunc_linprog
% Subfunctions: none
% MAT-files required: none
%
% See also: supportFunc_linprog

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   26-March-2026 (TL, delegate to supportFunc_linprog)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% simple check: empty polytope (fullspace)
if isempty(A) && isempty(Ae)
    switch type
        case 'upper'
            val = Inf; x = [];
        case 'lower'
            val = -Inf; x = [];
        case 'range'
            val = interval(-Inf,Inf); x = [];
    end
    return
end

% delegate to shared LP support function
[val,x] = supportFunc_linprog(dir,A,b,Ae,be,[],[],type);

% ------------------------------ END OF CODE ------------------------------
