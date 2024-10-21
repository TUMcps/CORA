function decodedArray = tree2evalstr(encodedArray,gp)
%TREE2EVALSTR Converts encoded tree expressions into math expressions that MATLAB can evaluate directly.
%
%   DECODEDARRAY = TREE2EVALSTR(ENCODEDARRAY,GP) processes the encoded
%   symbolic expressions in the cell array ENCODEDARRAY using the
%   information in the GP data structure and creates a cell array
%   DECODEDARRAY, each element of which is an evaluable Matlab string
%   expression.
%
%   Remarks: 
%
%   GPTIPS uses a encoded string representation of expressions which is not
%   very human readable, but is fairly compact and makes events like
%   crossover and mutation easier to handle. An example of such an
%   expression is a(b(x1,a(x4,x3)),c(x2,x1)). Before an expression can be
%   evaluated it is converted using TREE2EVALSTR to produce an evaluable
%   math expression, e.g. times(minus(x1,times(x4,x3)),plus(x2,x1)).
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
%
%   See also EVALFITNESS, TREEGEN, GPREFORMAT

%loop through active function list and replace with function names
for j=1:gp.nodes.functions.num_active
    encodedArray = strrep(encodedArray,gp.nodes.functions.afid(j),gp.nodes.functions.active_name_UC{j});
end

decodedArray = lower(encodedArray);