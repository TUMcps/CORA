function E = cartProd_(E,S,mode)
% cartProd_ - computes the Cartesian product of an ellipsoid and another
%    set or point
%
% Syntax:
%    E = cartProd_(E,S)
%    E = cartProd_(E,S,mode)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object, numeric
%    mode - type of approximation: 'exact', 'outer', 'inner'
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    E1 = ellipsoid([3 -1; -1 1],[1;0]);
%    E2 = ellipsoid([5 1; 1 2],[1;-1]);
%    E = cartProd(E1,E2); 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd

% Authors:       Victor Gassmann
% Written:       19-March-2021
% Last update:   02-June-2022 (VG, handle empty case)
%                17-April-2024 (TL, simplified ellipsoid-ellipsoid case)
% Last revision: 27-March-2023 (MW, rename cartProd_)
%                22-September-2024 (MW, restructure)

% ------------------------------ BEGIN CODE -------------------------------

switch mode
    case {'inner', 'exact'}
        throw(CORAerror('CORA:notSupported',...
            ['The function ''cartProd'' supports only type = ''outer'' ' ...
            ' for ellipsoid objects.']));

    case 'outer'
        if isa(E,'ellipsoid')
            if isa(S,'ellipsoid')
                % ellipsoid-ellipsoid case: Cartesian product of interval
                % over-approximations
                E = ellipsoid(cartProd_(interval(E),interval(S),'exact'));
                return
    
            elseif isnumeric(S) && iscolumn(S)
                % ellipsoid-numeric case
                E = ellipsoid(blkdiag(E.Q,zeros(length(S))), [E.q;S]);
                return
            end
    
        elseif isnumeric(E) && isa(S,'ellipsoid')
            % numeric-ellipsoid case
            E = ellipsoid(blkdiag(zeros(length(E)),S.Q), [E;S.q]);
            return

        end
end

% all other cases: throw error
throw(CORAerror('CORA:noops',E,S,mode));

end

% ------------------------------ END OF CODE ------------------------------
