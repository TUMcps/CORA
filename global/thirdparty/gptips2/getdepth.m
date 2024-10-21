function depth=getdepth(expr)
%GETDEPTH Returns the tree depth of an encoded tree expression.
%
%   DEPTH = GETDEPTH(EXPR) returns the DEPTH of the tree expression EXPR.
%
%   Remarks:
%
%   A single node is a tree of depth 1.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GETNUMNODES, GETCOMPLEXITY

%check to make sure expression is not actually a cell array of size 1
if iscell(expr)
    expr = expr{1};
end

%workaround for arity zero functions, which are represented by f()
%string replace '()' with the empty string
expr = strrep(expr,'()','');
open_br = strfind(expr,'(');
close_br = strfind(expr,')');
num_open = numel(open_br);

if ~num_open  %i.e. a single node
    depth = 1;
    return
elseif num_open == 1
    depth = 2;
    return
else
    %depth = max consecutive number of open brackets+1
    br_vec = zeros(1,numel(expr));
    br_vec(open_br) = 1;
    br_vec(close_br) = -1;
    depth = max(cumsum(br_vec)) + 1;
end