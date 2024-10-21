function comp = getcomplexity(expr)
%GETCOMPLEXITY Returns the expressional complexity of an encoded tree or a cell array of trees.
%
%   COMP = GETCOMPLEXITY(EXPR) finds the expressional complexity COMP of
%   the encoded GPTIPS expression EXPR or cell array of expressions EXPR.
%
%   Remarks:
%
%   "Expressional Complexity" is defined as the sum of nodes of all
%   sub-trees within a tree (i.e. as defined by Guido F. Smits and Mark
%   Kotanchek) in "Pareto-front exploitation in symbolic regression", pp.
%   283 - 299, Genetic Programming Theory and Practice II, 2004.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GETDEPTH, GETNUMNODES

if isa(expr,'char')
    
    comp = getcomp(expr);
    return;
    
elseif iscell(expr)
    
    numexpr = numel(expr);
    
    if numexpr < 1
        error('Cell array must contain at least one valid symbolic expression');
    else
        comp = 0;
        for i=1:numexpr
            comp = comp + getcomp(expr{i});
        end
    end
else
    error('Illegal argument');
end

function comp = getcomp(expr)
%GETCOMP Get complexity from a single tree.

%count leaf nodes (excluding zero arity functions).
%these have complexity measure of 1.
leafcount = numel(strfind(expr,'x')) + numel(strfind(expr,'['));

%iterate across remaining function
%nodes, getting the number of nodes of each subtree found
ind = strfind(expr,'(');
if isempty(ind)
    comp = leafcount;
    return;
end

ind = ind - 1;

comp = 0;
for i = 1:length(ind)
    [~,subtree] = extract(ind(i),expr);
    comp = comp + getnn(subtree);
end

comp = comp + leafcount;

function numnodes = getnn(expr)
%GETNN get number of nodes from a single tree

%number of nodes = number of open brackets + number of inputs + number of
%ERC constants
open_br = strfind(expr,'(');
open_sq_br = strfind(expr,'[');
inps = strfind(expr,'x');

num_open = numel(open_br);
num_const = numel(open_sq_br);
num_inps = numel(inps);
numnodes = num_open + num_const + num_inps;