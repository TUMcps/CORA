function [lsEq,lsIneq] = splitEquality(ls)
% splitEquality - split into level set with eqaulity constraints and level
%   set with inequality constraints
%
% Syntax:
%    [lsEq,lsIneq] = splitEquality(ls)
%
% Inputs:
%    ls - levelSet object
%
% Outputs:
%    lsEq - levelSet object containing only equality constraints
%    lsIneq - levelSet object containing only inequality constraints
%
% Example: 
%    syms x y
%    eq = [x^2 + y^2 - 1; x];
%    ls = levelSet(eq,[x;y],{'==';'<='});
%    [lsEq,lsIneq] = splitEquality(ls)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: levelSet

% Authors:       Niklas Kochdumper
% Written:       30-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    if ~iscell(ls.compOp)               % single comparison operator

        if strcmp(ls.compOp,'==')
            lsEq = ls; lsIneq = [];
        else
            lsIneq = ls; lsEq = [];
        end 

    else                                % multiple comparison operators

        indEq = find(strcmp(ls.compOp,'=='));
        indIneq  = setdiff(1:length(ls.compOp),indEq);

        if ~isempty(indEq)
            lsEq = levelSet(ls.eq(indEq),ls.vars,'==');
        else
            lsEq = [];
        end

        if ~isempty(indIneq)
            lsIneq = levelSet(ls.eq(indIneq),ls.vars,ls.compOp(indIneq));
        else
            lsIneq = [];
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
