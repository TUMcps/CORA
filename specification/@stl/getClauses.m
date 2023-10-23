function list = getClauses(obj,varargin)
% getClauses - list of clauses from STL formula in conjunctive normal form
%
% Syntax:
%    list = getClauses(obj)
%    list = getClauses(obj,type)
%
% Inputs:
%    obj - logic formula in conjunctive normal form (class stl)
%    type - conjuctive ('cnf') or disjunctive normal form('dnf')
%
% Outputs:
%    list - list of clauses
%
% Example: 
%    x = stl('x',2);
%    eq = ~(x(1) < 5 | x(2) < 3) | x(2) > 5;
%    eq_ = conjunctiveNormalForm(eq);
%    list = getClauses(eq_);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper, Benedikt Seidl
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    type = 'cnf';

    if nargin > 1
        type = varargin{1};
    end

    % get clauses 
    if ischar(type) && strcmp(type,'cnf')
        list = aux_getClausesCNF(obj);
    elseif ischar(type) && strcmp(type,'dnf')
        list = aux_getClausesDNF(obj);
    else
        throw(CORAerror('CORA:wrongValue','second',"'cnf' or 'dnf'"));
    end
end


% Auxiliary functions -----------------------------------------------------

function list = aux_getClausesCNF(obj)
% get clauses for conjunctive normal form

    if strcmp(obj.type,'&')
        list = [aux_getClausesCNF(obj.lhs);aux_getClausesCNF(obj.rhs)];
    else
        list = {obj};
    end
end

function list = aux_getClausesDNF(obj)
% get clauses for disjunctive normal form

    if strcmp(obj.type,'|')
        list = [aux_getClausesDNF(obj.lhs);aux_getClausesDNF(obj.rhs)];
    else
        list = {obj};
    end
end

% ------------------------------ END OF CODE ------------------------------
