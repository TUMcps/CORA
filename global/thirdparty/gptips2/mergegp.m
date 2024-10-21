function gp1 = mergegp(gp1,gp2,multirun)
%MERGEGP Merges two GP population structs into a new one.
%
%   GPNEW = MERGEGP(GP1,GP2) takes the results of 2 GPTIPS runs in GP1 and
%   GP2 and merges the populations into a new population GPNEW.
%
%   Remarks:
%
%   The merged populations do not have to have exactly the same run
%   settings in terms of population size, number of generations, selection
%   method type etc. but they do need to have the same function sets (i.e.
%   mathematical building blocks), inputs and be based on the same data set
%   (i.e. same training, validation and test sets). If you randomise the
%   allocation of training, test and validation sets in your config file,
%   then the results (in terms of fitness, R^2 etc.) will not be
%   meaningful.
%
%   Important:
%
%   The merged structure GPNEW can be used like any other GP results
%   structure (e.g. in POPBROWSER, RUNTREE etc.). However, many of the
%   fields in GPNEW will not necessarily apply to all the data in the new
%   structure. For instance, gp.runcontrol.pop_size may be different in GP2
%   to GP1 and gp.results.history and gp.state will certainly be different.
%   Note that GPNEW will always inherit these fields from GP1.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPMODELFILTER

disp('Merging GP populations ...');

if nargin < 3
   multirun = false; 
end

%do some tests to see if the GP structures are ostensibly OK to be merged.
if ~isequal(gp1.nodes.inputs,gp2.nodes.inputs)
    error('GP structures cannot be merged because the input nodes are not configured identically.');
end

if ~isequal(gp1.nodes.functions,gp2.nodes.functions)
    error('GP structures cannot be merged because the function nodes are not configured identically.');
end

%update the appropriate fields of gp1.
popsize1 = gp1.runcontrol.pop_size;
popsize2 = gp2.runcontrol.pop_size;
inds = popsize1+1:(popsize1+popsize2);

gp1.pop(inds) = gp2.pop(:);
gp1.fitness.returnvalues(inds) = gp2.fitness.returnvalues(:);
gp1.fitness.complexity(inds) = gp2.fitness.complexity(:);
gp1.fitness.nodecount(inds) = gp2.fitness.nodecount(:);
gp1.fitness.values(inds) = gp2.fitness.values(:);
gp1.runcontrol.pop_size = popsize1 + popsize2;

%check for existence of r2train r2val and r2test fields before merging
if isfield(gp1.fitness,'r2train') && isfield(gp2.fitness,'r2train')
    gp1.fitness.r2train(inds) = gp2.fitness.r2train(:);
end

if isfield(gp1.fitness,'r2val') && isfield(gp2.fitness,'r2val')
    gp1.fitness.r2val(inds) = gp2.fitness.r2val(:);
end

if isfield(gp1.fitness,'r2test') && isfield(gp2.fitness,'r2test')
    gp1.fitness.r2test(inds) = gp2.fitness.r2test(:);
end

%check to see if 'best' on validation data set is present
if isfield(gp1.results,'valbest') && isfield(gp2.results,'valbest')
    valPresent = true;
else
    valPresent = false;
end

%check to see if 'testbest' on test data set is present
if isfield(gp1.results,'testbest') && isfield(gp2.results,'testbest')
    testPresent = true;
else
    testPresent = false;
end

%find the best and valbest of gp1 and gp2 and get complexity values.
if gp1.fitness.complexityMeasure
    traincomp1 = gp1.results.best.complexity;
    traincomp2 = gp2.results.best.complexity;
    
    if valPresent
        valcomp1 = gp1.results.valbest.complexity;
        valcomp2 = gp2.results.valbest.complexity;
    end
    
    if testPresent
        testcomp1 = gp1.results.testbest.complexity;
        testcomp2 = gp2.results.testbest.complexity;
    end
    
else
    traincomp1 = gp1.results.best.nodecount;
    traincomp2 = gp2.results.best.nodecount;
    
    if valPresent
        valcomp1 = gp1.results.valbest.nodecount;
        valcomp2 = gp2.results.valbest.nodecount;
    end
    
    if testPresent
        testcomp1 = gp1.results.testbest.nodecount;
        testcomp2 = gp2.results.testbest.nodecount;
    end
end

if gp1.fitness.minimisation
    if gp2.results.best.fitness < gp1.results.best.fitness || ...
            ((gp2.results.best.fitness == gp1.results.best.fitness) && traincomp2 < traincomp1)
        gp1.results.best = gp2.results.best;
    end
    
    if valPresent
        if gp2.results.valbest.valfitness < gp1.results.valbest.valfitness || ...
                ((gp2.results.valbest.valfitness == gp1.results.valbest.valfitness) && valcomp2 < valcomp1)
            gp1.results.valbest = gp2.results.valbest;
        end
    end
    
    if testPresent
        if gp2.results.testbest.testfitness < gp1.results.testbest.testfitness || ...
                ((gp2.results.testbest.testfitness == gp1.results.testbest.testfitness) && testcomp2 < testcomp1)
            gp1.results.testbest = gp2.results.testbest;
        end
    end
    
else %maximisation
    if gp2.results.best.fitness > gp1.results.best.fitness || ...
            ((gp2.results.best.fitness == gp1.results.best.fitness) && traincomp2 < traincomp1)
        gp1.results.best = gp2.results.best;
    end
    
    if valPresent
        if gp2.results.valbest.valfitness > gp1.results.valbest.valfitness || ...
                ((gp2.results.valbest.valfitness == gp1.results.valbest.valfitness) && valcomp2 < valcomp1)
            gp1.results.valbest = gp2.results.valbest;
        end
    end
    
    if gp2.results.testbest.testfitness > gp1.results.testbest.testfitness || ...
            ((gp2.results.testbest.testfitness == gp1.results.testbest.testfitness) && testcomp2 < testcomp1)
        gp1.results.testbest = gp2.results.testbest;
    end
end

%merge info
if ~gp1.info.merged && ~gp2.info.merged
    gp1.info.mergedPopSizes = [popsize1 popsize2];
elseif gp1.info.merged && ~gp2.info.merged
    gp1.info.mergedPopSizes = [gp1.info.mergedPopSizes popsize2];
elseif ~gp1.info.merged && gp2.info.merged
    gp1.info.mergedPopSizes = [popsize1 gp2.info.mergedPopSizes];
elseif gp1.info.merged && gp2.info.merged
    gp1.info.mergedPopSizes = [gp1.info.mergedPopSizes gp2.info.mergedPopSizes];
end

%counter for number of times gp1 has experienced a merge
gp1.info.merged = gp1.info.merged + 1;
gp1.info.duplicatesRemoved = false;

%set stop time to time taken for all runs in multirun
if multirun
    gp1.info.stopTime = gp2.info.stopTime;
end