function gp = popbuild(gp)
%POPBUILD Build next population of individuals.
%
%   GP = POPBUILD(GP) uses the current population (stored in GP.POP) and
%   their fitnesses (stored in GP.FITNESS.VALUES) and, optionally,
%   complexities to create the next generation of individuals.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also INITBUILD

%initialise new population
newPop = cell(gp.runcontrol.pop_size,size(gp.userdata.ytrain,2));

%the number of new members to be constructed after elitism is accounted for
num2build = floor( (1 - gp.selection.elite_fraction) * gp.runcontrol.pop_size);
num2skim = gp.runcontrol.pop_size - num2build;

%parameter shortcuts
p_mutate = gp.operators.mutation.p_mutate;
p_direct = gp.operators.directrepro.p_direct;
maxDepth = gp.treedef.max_depth;
max_nodes = gp.treedef.max_nodes;
p_cross_hi = gp.genes.operators.p_cross_hi;
crossRate = gp.genes.operators.hi_cross_rate;
useMultiGene = gp.genes.multigene;
pmd = p_mutate + p_direct;
maxNodesInf = isinf(max_nodes);
maxGenes = gp.genes.max_genes;

%reset cache
if gp.runcontrol.usecache
    remove(gp.fitness.cache,gp.fitness.cache.keys);
end

%update gen counter
gp.state.count = gp.state.count + 1;

for i_y = 1:size(gp.userdata.ytrain,2)
    buildCount = 0;

    %loop until the required number of new individuals has been built.
    while buildCount < num2build;

        buildCount = buildCount + 1;

        %probabilistically select a genetic operator
        p_gen = rand;

        if p_gen < p_mutate  %select mutation
            eventType = 1;
        elseif p_gen < pmd   %direct reproduction
            eventType = 2;
        else                 %crossover
            eventType = 3;
        end


        %mutation (first select a individual, then a gene, then do ordinary mutate)
        if eventType == 1

            parentIndex = selection(gp);  %pick the population index of a parent individual using selection operator
            parent = gp.pop{parentIndex,i_y};

            if useMultiGene %if using multigene, extract a target gene
                numParentGenes = numel(parent);
                targetGeneIndex = ceil(rand * numParentGenes); %randomly select a gene from parent
                targetGene = parent{1,targetGeneIndex}; %extract it
            else
                targetGeneIndex = 1;
                targetGene = parent{1};
            end

            mutateSuccess = false;
            for loop = 1:10	%loop until a successful mutation occurs (max loops=10)

                mutatedGene = mutate(targetGene,gp);
                mutatedGeneDepth = getdepth(mutatedGene);

                if mutatedGeneDepth <= maxDepth

                    if maxNodesInf
                        mutateSuccess = true;
                        break;
                    end

                    mutatedGeneNodes = getnumnodes(mutatedGene);
                    if mutatedGeneNodes <= max_nodes
                        mutateSuccess = true;
                        break;
                    end
                end %end of constraint check

            end  %end of mutate for loop

            %if no success then use parent gene
            if ~mutateSuccess
                mutatedGene = targetGene;
            end

            %add the mutated individual to new pop
            parent{1,targetGeneIndex} = mutatedGene;
            newPop{buildCount,i_y} = parent;

            %direct reproduction operator
        elseif eventType == 2

            parentIndex = selection(gp);  %pick a parent
            parent = gp.pop{parentIndex,i_y};

            %copy to new population
            newPop{buildCount,i_y} = parent;

            %store fitness etc of copied individual if cache enabled
            if gp.runcontrol.usecache
                cachedData.complexity = gp.fitness.complexity(parentIndex,1);
                cachedData.returnvalues = gp.fitness.returnvalues{parentIndex,1};
                cachedData.value = gp.fitness.values(parentIndex,1);
                gp.fitness.cache(buildCount) = cachedData;
            end

            %crossover operator - can either pick 'high level' crossover
            %(crosses over entire genes with no tree alteration) or 'low level'
            % which crosses over individual genes at the tree level.
        elseif eventType == 3

            highLevelCross = false;

            if useMultiGene

                %select crossover type if multigene enabled
                if rand < p_cross_hi
                    highLevelCross = true;
                end

            end

            %Select Parents
            parentIndex = selection(gp);
            dad = gp.pop{parentIndex,i_y};
            numDadGenes = numel(dad);

            parentIndex = selection(gp);
            mum = gp.pop{parentIndex,i_y};
            numMumGenes = numel(mum);

            if highLevelCross
                if numMumGenes>1 || numDadGenes>1
                    %high level crossover (only use if either parent has more
                    %than one gene) This is modified/simplified from the
                    %version in GPTIPS v1. This version just chooses between 1
                    %and N genes to exchange randomly and independently for
                    %each parent with a probability set by
                    %gp.genes.operators.hi_cross_rate (where N =  number of
                    %genes in that individual). Like GPTIPS 1, this assumes
                    %that the order of genes within an individual is not
                    %important.

                    dadGeneSelectionInds = rand(1,numDadGenes) < crossRate;
                    mumGeneSelectionInds = rand(1,numMumGenes) < crossRate;

                    if ~any(dadGeneSelectionInds)
                        dadGeneSelectionInds(1,ceil(numDadGenes *rand)) = true;
                    end

                    if ~any(mumGeneSelectionInds)
                        mumGeneSelectionInds(1,ceil(numMumGenes *rand)) = true;
                    end

                    dadSelectedGenes = dad(dadGeneSelectionInds);
                    mumSelectedGenes = mum(mumGeneSelectionInds);

                    dadRemainingGenes = dad(~dadGeneSelectionInds);
                    mumRemainingGenes = mum(~mumGeneSelectionInds);

                    mumOffspring = [mumRemainingGenes dadSelectedGenes];
                    dadOffspring = [dadRemainingGenes mumSelectedGenes];

                    %curtail offspring longer than the max allowed number of genes
                    %before writing to new population (only write 1 if no space for 2 offspring)
                    newPop{buildCount,i_y} = mumOffspring(1:(min(end,maxGenes)));
                    buildCount = buildCount+1;

                    if buildCount <= num2build
                        newPop{buildCount,i_y} = dadOffspring(1:(min(end,maxGenes)));
                    end

                else
                    highLevelCross = false;
                end
            end

            %low level crossover: picks a random gene from each parent and
            %crosses them. The offspring replace the original target genes
            if ~highLevelCross

                if useMultiGene %if multigene then select target genes
                    dad_target_gene_num = ceil(rand*numDadGenes); %randomly select a gene from dad
                    mum_target_gene_num = ceil(rand*numMumGenes); %randomly select a gene from mum
                    dad_target_gene = dad{1,dad_target_gene_num};
                    mum_target_gene = mum{1,mum_target_gene_num};
                else
                    dad_target_gene_num = 1;
                    mum_target_gene_num = 1;
                    dad_target_gene = dad{1};
                    mum_target_gene = mum{1};
                end


                for loop = 1:10  %loop (max 10 times) until both children meet size constraints

                    %produce 2 offspring
                    [son,daughter] = crossover(mum_target_gene,dad_target_gene,gp);
                    son_depth = getdepth(son);

                    %check if both children meet size and depth constraints
                    %if true then break out of while loop and proceed
                    crossOverSuccess = false;
                    if son_depth <= maxDepth
                        daughter_depth = getdepth(daughter);
                        if daughter_depth <= maxDepth

                            if maxNodesInf
                                crossOverSuccess = true;
                                break;
                            end

                            son_nodes = getnumnodes(son);
                            if son_nodes <= max_nodes
                                daughter_nodes = getnumnodes(daughter);
                                if  daughter_nodes <= max_nodes
                                    crossOverSuccess = true;
                                    break;
                                end;
                            end
                        end
                    end

                end

                %if no success then re-insert parents
                if ~crossOverSuccess
                    son = dad_target_gene;
                    daughter = mum_target_gene;
                end

                %write offspring back to right gene positions in parents and write to population
                dad{1,dad_target_gene_num} = son;
                newPop{buildCount,i_y} = dad;

                buildCount = buildCount+1;

                if buildCount <= num2build
                    mum{1,mum_target_gene_num} = daughter;
                    newPop{buildCount,i_y} = mum;
                end

            end %end of if ~use_high

        end % end of op_type if

    end %end of ii-->num2build  for


    %skim off the existing elite individuals and stick them on the end of the
    %new population. When sorting there may be many individuals with the same
    %fitness. However, some of these have fewer nodes/lower complexity than
    %others, so skim off the first one as the existing member with the fewest
    %nodes. This should exert a degree of parsimony pressure.

    %get indices of best num2skim individuals
    [~,sortIndex] = sort(gp.fitness.values);

    %if maximising need to flip vectors
    if ~gp.fitness.minimisation
        sortIndex = flipud(sortIndex);
    end

    for g=1:num2skim

        oldIndex = sortIndex(g);

        %for the first individual to skim off, pick the 'best' member with the
        %lowest complexity
        if g==1
            bestInds = find( gp.fitness.values == gp.fitness.values(oldIndex) );
            [~,oldIndex] = sort(gp.fitness.complexity(bestInds));
            oldIndex = bestInds(oldIndex(1)); %if more than one with same complexity just pick the first
        end

        newIndex = num2build + g;
        copiedIndividual = gp.pop{oldIndex,i_y};

        %cache fitness, complexity etc.
        if gp.runcontrol.usecache
            cachedData.complexity = gp.fitness.complexity(oldIndex,1);
            cachedData.returnvalues = gp.fitness.returnvalues{oldIndex,1};
            cachedData.value = gp.fitness.values(oldIndex,1);
            gp.fitness.cache(newIndex) = cachedData;
        end

        newPop{newIndex,i_y} = copiedIndividual;

    end
end

gp.pop = newPop;