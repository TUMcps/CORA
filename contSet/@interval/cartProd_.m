function res = cartProd_(I,S,varargin)
% cartProd_ - returns the Cartesian product of two intervals
%
% Syntax:
%    res = cartProd_(I,S)
%
% Inputs:
%    I - interval object
%    S - interval object
%
% Outputs:
%    res - Cartesian product of intervals
%
% Example: 
%    I1 = interval([-1;-3],[1;6]);
%    I2 = interval([-2;1], [4;2]);
%    res = cartProd(I1,I2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd, zonotope/cartProd_

% Authors:       Mark Wetzlinger
% Written:       18-September-2019
% Last update:   24-September-2019
%                05-May-2020 (standardized error message)
% Last revision: 27-March-2023 (MW, rename cartProd_)

% ------------------------------ BEGIN CODE -------------------------------

% first or second set is interval
if isa(I,'interval')

    % different cases for different set representations
    if isa(S,'interval')

        if representsa_(I,'emptySet',eps)
            res = S;
        elseif representsa_(S,'emptySet',eps)
            res = I;
        elseif size(I,2) == 1 && size(S,2) == 1 % vertcat
            res = interval([I.inf;S.inf],[I.sup;S.sup]);
        elseif size(I,1) == 1 && size(S,1) == 1 % horzcat
            res = interval([I.inf,S.inf],[I.sup,S.sup]);
        else
            throw(CORAerror('CORA:dimensionMismatch',I,S));
        end

    elseif isnumeric(S)

        res = interval([I.inf;S],[I.sup;S]);

    % different cases for different set representations
    elseif isa(S,'zonotope') 
        res = cartProd_(zonotope(I),S,'exact');
    elseif isa(S,'conZonotope')
        res = cartProd_(conZonotope(I),S,'exact');
    elseif isa(S,'zonoBundle')
        res = cartProd(zonoBundle(I),S,'exact');
    elseif isa(S,'polytope')
        res = cartProd_(polytope(I),S,'exact');
    elseif isa(S,'polyZonotope')
        res = cartProd_(polyZonotope(I),S,'exact');
    elseif isa(S,'conPolyZono')
        res = cartProd_(conPolyZono(I),S,'exact');
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',I,S));
    end

else

    if isnumeric(I)
        res = interval([I;S.inf],[I;S.sup]);
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',I,S));
    end  
end

% ------------------------------ END OF CODE ------------------------------
