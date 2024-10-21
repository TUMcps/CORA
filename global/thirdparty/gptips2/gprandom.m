function gp = gprandom(gp,seed)
%GPRANDOM Sets random number generator seed according to system clock or user seed.
%
%   (c) Dominic Searson 2009-2015
%
%   GPTIPS 2

if nargin < 2 || isempty(seed)
    seed = sum(100*clock);
end

gp.info.PRNGseed = seed;

if verLessThan('matlab', '7.7.0');
    rand('twister', gp.info.PRNGseed);
elseif verLessThan('matlab','8.1')
    s = RandStream.create('mt19937ar','seed',gp.info.PRNGseed);
    RandStream.setDefaultStream(s);
else
    s = RandStream.create('mt19937ar','seed',gp.info.PRNGseed);
    RandStream.setGlobalStream(s);
end