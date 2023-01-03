function res = cartProd(I,S)
% cartProd - returns the Cartesian product of two intervals
%
% Syntax:  
%    res = cartProd(I,S)
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
% See also: zonotope/cartProd

% Author:       Mark Wetzlinger
% Written:      18-Sep-2019
% Last update:  24-Sep-2019
%               05-May-2020 (standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[resFound,vars] = pre_cartProd('interval',I,S);

% check premature exit
if resFound
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    % assign variables
    I = vars{1}; S = vars{2};
end


% first or second set is interval
if isa(I,'interval')

    % different cases for different set representations
    if isa(S,'interval')

        if isempty(I)
            res = S;
        elseif isempty(S)
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
        res = cartProd(zonotope(I),S);
    elseif isa(S,'conZonotope')
        res = cartProd(conZonotope(I),S);
    elseif isa(S,'zonoBundle')
        res = cartProd(zonoBundle(I),S);
    elseif isa(S,'mptPolytope')
        res = cartProd(mptPolytope(I),S);
    elseif isa(S,'polyZonotope')
        res = cartProd(polyZonotope(I),S);
    elseif isa(S,'conPolyZono')
        res = cartProd(conPolyZono(I),S);
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

%------------- END OF CODE --------------