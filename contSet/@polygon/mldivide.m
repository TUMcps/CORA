function pgon_out = mldivide(pgon,S)
% mldivide - overloads the '\' operator
%
% Syntax:
%    pgon_out = mldivide(pgon, S)
%    pgon_out = pgon \ S
%
% Inputs:
%    pgon - polygon
%    S - contSet
%
% Outputs:
%    pgon_out - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
equalDimCheck(pgon,S);
inputArgsCheck({ ...
    {pgon,'att','polygon'}, ...
    {S,'att','contSet'}, ...
});

% convert S to polygon
S_pgon = polygon(S);

% use 'subtract' of polyshape
pgon_out = polygon(subtract(pgon.set,S_pgon.set));

end

% ------------------------------ END OF CODE ------------------------------
