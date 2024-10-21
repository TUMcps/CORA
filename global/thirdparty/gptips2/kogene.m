function treestrs = kogene(treestrs, knockout)
%KOGENE Knock out genes from a cell array of tree expressions.
%
%   TREESTRS = KOGENE(TREESTRS, KNOCKOUT) removes genes from the cell array
%   TREESTRS with a boolean vector KNOCKOUT which must be the same length
%   as TREESTRS. GENES corresponding to a an entry of 'true' are removed.
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2

if numel(knockout) ~= numel(treestrs)
    error('Knockout must be a vector the same length as the number of genes in the supplied individual');
end

knockout = logical(knockout);

treestrs(knockout) = [];

if isempty(treestrs)
    error('Cannot perform knockout. There are no genes left in the selected individual.');
end
