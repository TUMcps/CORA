function [symbolicResult, parallelResult, statsResult] = gptoolboxcheck
%GPTOOLBOXCHECK Checks if certain toolboxes are installed and licensed.
%
%   Checks for the Stats, Parallel Computing and Symbolic Math toolboxes.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPINIT

v = ver;
[toolboxes{1:length(v)}] = deal(v.Name);
symbolicResult = ismember('Symbolic Math Toolbox',toolboxes) && license('test','symbolic_toolbox');
parallelResult = ismember('Parallel Computing Toolbox',toolboxes) && license('test','distrib_computing_toolbox');
statsResult = ( ismember('Statistics Toolbox',toolboxes) || ismember('Statistics and Machine Learning Toolbox',toolboxes) ) ...
    && license('test','statistics_toolbox'); 

