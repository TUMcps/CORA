function gp=initbuild(gp)
%INITBUILD Generate an initial population of GP individuals.
%
%   GP = INITBUILD(GP) creates an initial population using the parameters
%   in the structure GP. Various changes to fields of GP are made.
%
%   Remarks:
%   
%   Each individual in the population can contain 1 or more separate trees
%   (if more than 1 then each tree is referred to as a gene of the
%   individual). This is set by the user in the config file in the field
%   GP.GENES.MULTIGENE.
%
%   Each individual is a cell array. You can use cell addressing to
%   retrieve the individual genes. E.g. GP.POP{3}{2} will return the 2nd
%   gene in the third individual.
%
%   Trees are (by default) constructed using a probabilistic version of
%   Koza's ramped 1/2 and 1/2 method. E.g. if maximum tree depth is 5 and
%   population is 100 then, on average, 20 will be generated at each depth
%   (1/2 using 'grow' and 1/2 using 'full').
%
%   Multiple copies of genes in an individual are disallowed when
%   individuals are created. There is, however, no such restriction imposed
%   on future generations.
%
%   (c) Dominic Searson 2009-2015
%
%   GPTIPS 2
%
%   See also POPBUILD, TREEGEN

% Extract temp variables from the gp structure
popSize = gp.runcontrol.pop_size;
maxNodes = gp.treedef.max_nodes;
maxGenes = gp.genes.max_genes;

%override any gene settings if using single gene gp
if ~gp.genes.multigene
    maxGenes = 1;
end

%initialise vars
gp.pop = cell(popSize,size(gp.userdata.ytrain,2));
for i_y = 1: size(gp.userdata.ytrain,2)
    numGenes = 1;

    %building process
    for i=1:popSize %loop through population

        %randomly pick num of genes in individual
        if maxGenes > 1
            numGenes = ceil(rand*maxGenes);
        end

        individ = cell(1,(numGenes)); %construct empty individual

        for z = 1:numGenes %loop through genes in each individual and generate a tree for each

            %generate gene z and check that genes 1...z-1 are different
            geneLoopCounter = 0;
            while true

                geneLoopCounter = geneLoopCounter + 1;

                %generate a trial tree for gene z
                temp = treegen(gp);
                numnodes = getnumnodes(temp);

                if numnodes <= maxNodes

                    copyDetected = false;

                    if z > 1 %check previous genes for copies

                        for j = 1:z-1;
                            if strcmp(temp,individ{1,j})
                                copyDetected = true;
                                break
                            end
                        end

                    end

                    if ~copyDetected
                        break
                    end

                end %max nodes check

                %display a warning if having difficulty building trees due to
                %constraints
                if ~gp.runcontrol.quiet && geneLoopCounter > 10
                    disp('initbuild: iterating tree build loop because of uniqueness constraints.');
                end

            end %while loop

            individ{1,z} = temp;

        end %gene loop

        %write new individual to population cell array
        gp.pop{i,i_y} = individ;
    end
end
