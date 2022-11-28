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
%    x = stl('x',2)
%    eq = ~(x(1) < 5 | x(2) < 3) | x(2) > 5;
%    eq_ = conjunctiveNormalForm(eq);
%    list = getClauses(eq_);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper
% Written:      9-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'cnf';

    if nargin > 1
        type = varargin{1};
    end

    % get clauses 
    if ischar(type) && strcmp(type,'cnf')
        list = getClausesCNF(obj);
    elseif ischar(type) && strcmp(type,'dnf')
        list = getClausesDNF(obj);
    else
        throw(CORAerror('CORA:wrongValue','second',"'cnf' or 'dnf'"));
    end
end


% Auxiliary Functions -----------------------------------------------------

function list = getClausesCNF(obj)
% get clauses for conjunctive normal form

    if strcmp(obj.type,'&')
        list = [getClausesCNF(obj.lhs);getClausesCNF(obj.rhs)];
    else
        list = {obj};
    end
end

function list = getClausesDNF(obj)
% get clauses for disjunctive normal form

    if strcmp(obj.type,'|')
        list = [getClausesDNF(obj.lhs);getClausesDNF(obj.rhs)];
    else
        list = {obj};
    end
end

%------------- END OF CODE --------------