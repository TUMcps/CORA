function gp = evalfitness_par(gp)
%EVALFITNESS_PAR Calls the user specified fitness function (parallel version).
%
%   GP = EVALFITNESS_PAR(GP) evaluates the the fitnesses of individuals
%   stored in the GP structure and updates various other fields of GP
%   accordingly.
%
%   This has the same functionality as EVALFITNESS but makes use of the
%   Mathworks Parallel Computing Toolbox to distribute the fitness
%   computations across multiple cores.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also TREE2EVALSTR, EVALFITNESS

popSize = gp.runcontrol.pop_size;
evalstrs = cell(popSize,1);
complexityMeasure = gp.fitness.complexityMeasure;
complexities = zeros(popSize,1);
fitfun = gp.fitness.fitfun;
fitvals = zeros(popSize,1);
returnvals = cell(popSize,1);

if gp.runcontrol.usecache
    usecache = false; %true;
else
    usecache = false;
end

parfor i = 1:popSize;
    
    %assign copy of gp inside parfor loop to a temp struct
    tempgp = gp;
    
    %update state of temp variable to index of the individual that is about to
    %be evaluated
    tempgp.state.current_individual = i;
    
    if usecache && tempgp.fitness.cache.isKey(i)
        
        cache = tempgp.fitness.cache(i);
        complexities(i) = cache.complexity;
        fitvals(i) = cache.value;
        returnvals{i} = cache.returnvalues;
        
    else

        evalstr_all = cell(tempgp.genes.max_genes, size(tempgp.userdata.ytrain,2));
        [evalstr_all{:}] = deal('0');
        compl = 0;
        for i_y = 1:size(gp.userdata.ytrain,2)
            evalstr = tree2evalstr(tempgp.pop{i,i_y},gp);
            evalstr_all(1:length(evalstr),i_y) = evalstr;

            %store complexity of individual (either number of nodes or tree
            %expressional complexity)
            if complexityMeasure
                compl = compl + getcomplexity(tempgp.pop{i,i_y});
            else
                compl = compl + getnumnodes(tempgp.pop{i,i_y});
            end
        end
        complexities(i) = compl;

        %process coded trees into evaluable matlab expressions
        evalstrs{i} = evalstr_all;
        s_id = string(i);

        %evaluate gp individual using fitness function
        [fitness,tempgp] = feval(fitfun,evalstrs{i},tempgp);
        returnvals{i} = tempgp.fitness.returnvalues{i};
        fitvals(i) = fitness;
    end
end %end of parfor loop

%attach returned values to original GP structure
gp.fitness.values = fitvals;
gp.fitness.returnvalues = returnvals;
gp.fitness.complexity = complexities;
gp.state.current_individual = popSize;