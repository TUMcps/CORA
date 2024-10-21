function gp = evalfitness(gp)
%EVALFITNESS Calls the user specified fitness function.
%
%   GP = EVALFITNESS(GP) evaluates the the fitnesses of individuals stored
%   in the GP structure and updates various other fields of GP accordingly.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also TREE2EVALSTR, EVALFITNESS_PAR

%check parallel mode.
if gp.runcontrol.parallel.enable && gp.runcontrol.parallel.ok
    gp = evalfitness_par(gp);
    return;

    %regular version
else
    for i = 1:size(gp.pop,1)%gp.runcontrol.pop_size

        gp.state.current_individual = i;

        %retrieve values if cached
        if gp.runcontrol.usecache && gp.fitness.cache.isKey(i)
            cache = gp.fitness.cache(i);
            gp.fitness.complexity(i) = cache.complexity;
            gp.fitness.values(i) = cache.value;
            gp.fitness.returnvalues{i} = cache.returnvalues;

        else
            %preprocess cell array of string expressions into a form that
            %Matlab can evaluate
            evalstr_all = cell(gp.genes.max_genes, size(gp.userdata.ytrain,2));
            [evalstr_all{:}] = deal('0');
            compl = 0;
            for i_y = 1:size(gp.userdata.ytrain,2)
                evalstr = tree2evalstr(gp.pop{i,i_y},gp);
                evalstr_all(1:length(evalstr),i_y) = evalstr;

                %store complexity of individual (either number of nodes or tree
                %expressional complexity)
                if gp.fitness.complexityMeasure
                    compl = compl + getcomplexity(gp.pop{i, i_y});
                else
                    compl = compl + getnumnodes(gp.pop{i,i_y});
                end
            end
            gp.fitness.complexity(i) = compl;

            if nargin(gp.fitness.fitfun)==3
                [fitness,gp] = feval(gp.fitness.fitfun,evalstr_all,gp,i);
            else
                [fitness,gp] = feval(gp.fitness.fitfun,evalstr_all,gp);
            end
            gp.fitness.values(i) = fitness;
        end
    end
end
end