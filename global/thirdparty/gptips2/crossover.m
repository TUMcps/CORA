function [son,daughter] = crossover(mum,dad,gp)
%CROSSOVER Sub-tree crossover of encoded tree expressions to produce 2 new ones.
%
%   [SON,DAUGHTER] = CROSSOVER(MUM,DAD,GP) uses standard subtree crossover
%   on the expressions MUM and DAD to produce the offspring expressions SON
%   and DAUGHTER.
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
%
%   See also MUTATE

%select random crossover nodes in mum and dad expressions
m_position = picknode(mum,0,gp);
d_position = picknode(dad,0,gp);

%extract main and subtree expressions
[m_main,m_sub] = extract(m_position,mum);
[d_main,d_sub] = extract(d_position,dad);

%combine to form 2 new GP trees
daughter = strrep(m_main,'$',d_sub);
son = strrep(d_main,'$',m_sub);