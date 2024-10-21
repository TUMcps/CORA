function gp = gpterminate(gp)
%GPTERMINATE Check for early termination of run.
%
%   GP = GPTERMINATE(GP) checks if the GPTIPS run should be terminated at
%   the end of the current generation according to a fitness criterion or
%   run timeout.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPFINALISE, GPTIC, GPTOC

%check if fitness termination criterion met
if gp.fitness.terminate
    
    if gp.fitness.minimisation
        if gp.results.best.fitness <= gp.fitness.terminate_value
            gp.state.terminate = true;
            gp.state.terminationReason = ['Fitness criterion acheived <= ' ...
                num2str(gp.fitness.terminate_value)];
        end
    else
        if gp.results.best.fitness >= gp.fitness.terminate_value;
            gp.state.terminate = true;
            gp.state.terminationReason = ['Fitness criterion acheived >= ' ...
                num2str(gp.fitness.terminate_value)];
        end
    end
end

%run timeout
if ~gp.state.terminate && (gp.state.runTimeElapsed  >= gp.runcontrol.timeout)
    gp.state.terminate = true;
    gp.state.terminationReason = ['Timeout >= ' num2str(gp.runcontrol.timeout) ' sec'];
end