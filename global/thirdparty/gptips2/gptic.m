function gp = gptic(gp)
%GPTIC Updates the running time of this run.
%
%   GP = GPTIC(GP) begins the start timer for the current generation. 
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPTOC

gp.state.tic = tic;