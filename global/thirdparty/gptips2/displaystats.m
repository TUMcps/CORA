function displaystats(gp)
%DISPLAYSTATS Displays run stats periodically.
%
%   DISPLAYSTATS(GP) updates the screen with run stats at the interval
%   specified in GP.RUN.VERBOSE
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
%
%   See also UPDATESTATS

%only display info if required
if  ~gp.runcontrol.verbose || gp.runcontrol.quiet || ...
        mod(gp.state.count-1,gp.runcontrol.verbose)
    return
end

gen = gp.state.count - 1;

%multiple independent run count display
if gp.runcontrol.runs > 1  && gen == 0
   disp(['Run ' int2str(gp.state.run) ' of ' int2str(gp.runcontrol.runs)]);
end

disp(['Generation ' num2str(gen)]);
disp(['Best fitness:    ' num2str(gp.results.best.fitness)]);
disp(['Mean fitness:    ' num2str(gp.state.meanfitness)]);

if gp.fitness.complexityMeasure
    disp(['Best complexity: ' num2str(gp.results.best.complexity)]);
else
    disp(['Best node count: ' num2str(gp.results.best.nodecount)]);
end

%display inputs in "best training" individual, if enabled.
if gp.runcontrol.showBestInputs
    hitvec = gpmodelvars(gp,'best');
    inputs = num2str(find(hitvec));
    inputs = strtrim(cellstr(inputs));
    pat = '(\d+)';
    inputs = regexprep(inputs,pat,'x$1 '); %note: trailing space essential
    inputs = [inputs{:}]; % converts char cell array to a single string
    disp(['Inputs in best individual: ' inputs]);
end

%display inputs in "best validation" individual, if enabled.
if gp.runcontrol.showValBestInputs && isfield(gp.userdata,'xval') ...
        && ~isempty(gp.userdata.xval) && isfield(gp.results,'valbest')
    hitvec = gpmodelvars(gp,'valbest');
    inputs = num2str(find(hitvec));
    inputs = strtrim(cellstr(inputs));
    pat='(\d+)';
    inputs = regexprep(inputs,pat,'x$1 '); %note: trailing space essential
    inputs = [inputs{:}]; % converts char cell array to a single string
    disp(['Inputs in best validation individual: ' inputs]);
end

disp(' ');