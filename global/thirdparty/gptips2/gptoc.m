function gp = gptoc(gp)
%GPTOC Updates the running time of this run.
%
%   GP = GPTOC(GP) updates the running time for the current run. 
%
%   Copyright (c) 2014 Dominic Searson 
%
%   GPTIPS 2
%
%   See also GPTIC

gp.state.runTimeElapsed = gp.state.runTimeElapsed + toc(gp.state.tic);

