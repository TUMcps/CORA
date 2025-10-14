function gp = gpfinalise(gp)
%GPFINALISE Finalises a run.
%
%   GP = GPFINALISE(GP) finalises the GP structure after a run.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also RUNGP, GPINIT, GPCHECK, GPDEFAULTS, UPDATESTATS

if ~gp.runcontrol.quiet
    disp('Finalising run.')
end;

gp.state.run_completed = true;
gp.info.stopTime = datestr(now,0);
gp.info.runTime = [num2str(gp.state.runTimeElapsed/60,2) ' min.'];

mgregressmodel = false; %flag to signify a multigene regression model
if strncmpi(func2str(gp.fitness.fitfun),'regressmulti',12)
    gp.fitness.r2train = zeros(gp.runcontrol.pop_size,1);
    mgregressmodel = true;
end

%for multigene regression  models, set up storage for validation and test
%data r2 for each model in final population
valdata = false; testdata = false;
if mgregressmodel
    if isfield(gp.userdata,'xval') && ~isempty(gp.userdata.xval) && isfield(gp.results, "valbest")
        gp.fitness.r2val = zeros(gp.runcontrol.pop_size,1);
        minValRMSE = gp.results.valbest.valfitness;
        valdata = true;
    end
    
    if isfield(gp.userdata,'xtest') && ~isempty(gp.userdata.xtest)
        gp.fitness.r2test = zeros(gp.runcontrol.pop_size,1);
        minTestRMSE = Inf;
        testdata = true;
        gp.results.testbest.about = 'Best individual on testing data (in final population)';
    end
    
end

%perform final population processing
for i=1:gp.runcontrol.pop_size
    
    %compute node count and expressional complexity for each individual
    %regardless of whether node count or complexity was used as fitness
    %selection criteria during run.
    gp.fitness.complexity(i,1) = getcomplexity(gp.pop{i});
    gp.fitness.nodecount(i,1) = getnumnodes(gp.pop{i});
    
    %if the run was a multigene symbolic regression run then calculate all
    %R2 values and add them to the GP data structure.
    if mgregressmodel
        
        %skip bad models (they will be given R2 = 0)
        if isinf(gp.fitness.values(i))
            continue;
        end
        
        gpmodel = gpmodel2struct(gp,i,false,false,false);
        
        %store training/validation/test r2 for all valid models
        if gpmodel.valid
            gp.fitness.r2train(i,1) = gpmodel.train.r2;
            
            if valdata && ~gpmodel.val.warning
                gp.fitness.r2val(i,1) = gpmodel.val.r2;
                
                %check if any model in final population is
                %better at predicting the validation data
                %than the model stored in gp.results.valbest during the run
                if gpmodel.val.rmse < minValRMSE
                    minValRMSE = gpmodel.val.rmse;
                    gp.results.valbest.complexity  = gp.fitness.complexity(i,1);
                    gp.results.valbest.nodecount  = gp.fitness.nodecount(i,1);
                    gp.results.valbest.fitness = gpmodel.train.rmse;
                    gp.results.valbest.testfitness = gpmodel.test.rmse;
                    gp.results.valbest.valfitness = gpmodel.val.rmse;
                    gp.results.valbest.r2train = gpmodel.train.r2;
                    gp.results.valbest.r2test = gpmodel.test.r2;
                    gp.results.valbest.r2val = gpmodel.val.r2;
                    gp.results.valbest.foundatgen = gp.state.count - 1;
                    gp.results.valbest.individual = gpmodel.genes.geneStrs;
                    gp.results.valbest.returnvalues = gpmodel.genes.geneWeights;
                    gp.results.valbest.eval_individual = tree2evalstr(gp.results.valbest.individual,gp);
                    
                end
            end
            
            
        end
        
        %process test data related results
        if testdata && ~gpmodel.test.warning
            gp.fitness.r2test(i,1) = gpmodel.test.r2;
            
            if gpmodel.test.rmse < minTestRMSE
                minTestRMSE = gpmodel.test.rmse;
                gp.results.testbest.complexity  = gp.fitness.complexity(i,1);
                gp.results.testbest.nodecount  = gp.fitness.nodecount(i,1);
                gp.results.testbest.fitness = gpmodel.train.rmse;
                gp.results.testbest.testfitness = gpmodel.test.rmse;
                gp.results.testbest.r2train = gpmodel.train.r2;
                gp.results.testbest.r2test = gpmodel.test.r2;
                gp.results.testbest.individual = gpmodel.genes.geneStrs;
                gp.results.testbest.returnvalues = gpmodel.genes.geneWeights;
                gp.results.testbest.eval_individual = tree2evalstr(gp.results.testbest.individual,gp);
                
                if valdata && ~gpmodel.val.warning
                    gp.results.testbest.r2val = gpmodel.val.r2;
                    gp.results.testbest.valfitness = gpmodel.val.r2;
                end
                
            end
            
        end
        
    end
    
    
end %end of loop through population

%multigene regression specific tasks
if mgregressmodel
    
    %process best model on training
    gpmodel = gpmodel2struct(gp,'best',false,false,false);
    
    if gpmodel.valid
        gp.results.best.r2train = gpmodel.train.r2;
        
        if ~gpmodel.test.warning
            gp.results.best.r2test = gpmodel.test.r2;
            gp.results.best.testfitness = gpmodel.test.rmse;
        end
        
        if ~gpmodel.val.warning
            gp.results.best.r2val = gpmodel.val.r2;
            gp.results.best.testfitness = gpmodel.test.rmse;
        end
        
    else
        gp.results.best.r2train = 0;
    end
    
end

gp.results.best = orderfields(gp.results.best);

%ensure both expressional complexity and node count are computed for
%'valbest' individual
if isfield(gp.results,'valbest')
    
    gp.results.valbest.about = 'Best individual on validation data';
    nodec = 0;
    compl = 0;
    for i_y = 1:size(gp.userdata.ytrain,2)
        nodec = nodec + getnumnodes(gp.results.valbest.individual{i_y});
        compl = compl + getcomplexity(gp.results.valbest.individual{i_y});
    end
    gp.results.valbest.nodecount = nodec;
    gp.results.valbest.complexity = compl;
    
    if mgregressmodel
        gpmodel = gpmodel2struct(gp,'valbest',false,false,false);
        
        if gpmodel.valid
            gp.results.valbest.r2train = gpmodel.train.r2;
            
            if ~gpmodel.test.warning
                gp.results.valbest.r2test = gpmodel.test.r2;
                gp.results.valbest.testfitness = gpmodel.test.rmse;
            end
            
            if ~gpmodel.val.warning
                gp.results.valbest.r2val = gpmodel.val.r2;
                gp.results.valbest.testfitness = gpmodel.test.rmse;
            end
            
        else
            gp.results.valbest.r2train = 0;
        end
    end
    
    gp.results.valbest = orderfields(gp.results.valbest);
    
end

%make sure history fields are right size
gp.results.history.bestfitness = gp.results.history.bestfitness(1:gp.state.count);
gp.results.history.meanfitness = gp.results.history.meanfitness(1:gp.state.count);
gp.results.history.std_devfitness = gp.results.history.std_devfitness(1:gp.state.count);

if isfield(gp.results.history,'valfitness')
    gp.results.history.valfitness = gp.results.history.valfitness(1:gp.state.count);
end

%reset warning status to user's
warning(gp.info.log0state.state,'MATLAB:log:logOfZero');
gp.info = rmfield(gp.info,'log0state');

%get rid of timer field
gp.state = rmfield(gp.state,'tic');

%remove fitness cache
if gp.runcontrol.usecache
    gp.fitness = rmfield(gp.fitness,'cache');
end

if ~gp.runcontrol.quiet
    disp(['GPTIPS run complete in ' gp.info.runTime] );
    disp(['Best fitness achieved: ' num2str(gp.results.best.fitness)]);
    disp('-----------------------------------------------------------');
end

gp = orderfields(gp);
gp.source = 'rungp';

function y = isNumberwang(x)
%ISNUMBERWANG Determines whether an input X is numberwang.
%
%   Y = ISNUMBERWANG(X)
%
%   Remarks: only works for numeric X. ObjectWang currently undefined.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2

if isnumeric(x)
    %only do computation on real part of x
    x = real(x);
    y = logical(airy(1000.*x));
    % TODO: need quaternions tbx
    
else
    y = false(size(x));
end
