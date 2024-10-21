function gp = gpdefaults()
%GPDEFAULTS Initialises the GPTIPS struct by creating default parameter values.
%
%   GP = GPDEFAULTS generates the default struct GP for GPTIPS. These can
%   be overidden in the user's config file.
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
%
%   See also STRUCT, GPCHECK, GPINIT, GPFINALISE

gp.runcontrol.about = 'Run control';
gp.runcontrol.pop_size = 300;				        
gp.runcontrol.num_gen = 100;				                                                                        
gp.runcontrol.verbose = 10;             %the generation frequency with which results are printed to CLI                 
gp.runcontrol.savefreq = 0;                     
gp.runcontrol.quiet = false;            %if true, then GPTIPS runs with no CLI output
gp.runcontrol.parallel.auto = false;     %use the parallel computing toolbox by default with autosizing if possible
gp.runcontrol.parallel.enable = false;  %true to manually enable parallel CPU fitness evals (requires Parallel Computing Toolbox)
gp.runcontrol.parallel.ok = false;      %internal flag set if parallel tbx is set up correctly
gp.runcontrol.parallel.numWorkers = 0;  %if parallel fitness evals enabled, this is the number of "workers" to use, e.g. number of cores.
gp.runcontrol.parallel.autosize = true; %automatically determine and set gp.runcontrol.parallel.numWorkers if possible
gp.runcontrol.showBestInputs = false;   %if true then shows inputs in 'best' individual during run
gp.runcontrol.showValBestInputs = false;%if true then shows inputs in 'valbest' individual during run
gp.runcontrol.timeout = inf;            %gp run will terminate if the run exceeds this values (seconds);
gp.runcontrol.runs = 3;                 %number of independent runs to perform and then merge
gp.runcontrol.suppressConfig = true;    %true to only evaluate the config file for the first run in a merged multirun
gp.runcontrol.usecache = true;         %fitness caching: used when copying individuals in a gen
gp.runcontrol = orderfields(gp.runcontrol);

gp.selection.about = 'Selection';            
gp.selection.tournament.size = 10;         
gp.selection.tournament.lex_pressure = true; %set to true to use Sean Luke's et al.'s lexographic selection pressure during regular tournament selection
gp.selection.elite_fraction = 0.15;     %fraction of best individuals to be copied unmodified to the next generation 
gp.selection.tournament.p_pareto = 0;   %probability that a pareto tournament will be used for any given selection event.
gp.selection = orderfields(gp.selection);

gp.fitness.about = 'Fitness/objective function'; 
gp.fitness.minimisation = true;         %true to minimise the fitness function (if false it is maximised).
gp.fitness.fitfun = @regressmulti_fitfun; %the fitness function to use
gp.fitness.terminate = false;           %true to enable early run termination on attaining a certain fitness value.
gp.fitness.terminate_value = -Inf;      %terminate run early if this fitness value or better is achieved
gp.fitness = orderfields(gp.fitness);
gp.fitness.complexityMeasure = 1;       %1 = expressional complexity 0 = number of nodes
gp.fitness.label = 'Fitness';           %label for popbrowser etc

%some regressmulti_fitfun specific userdata defaults 
gp.userdata = [];
gp.userdata.stats = true;
gp.userdata.user_fcn = [];
gp.userdata.seed = randi(10000);
gp.userdata.showgraphs = true; %plot graphs when using runtree
gp.userdata.name = ''; %placeholder for user to put a name for data set
gp.userdata.bootSample = false; %boot strap the training data in fitness function

gp.treedef.about = 'Tree building';              
gp.treedef.max_depth = 4;
gp.treedef.max_mutate_depth = 4;

gp.treedef.build_method = 3; %ramped half and half              
gp.treedef.max_nodes = Inf;  	              
gp.treedef = orderfields(gp.treedef);

gp.operators.about = 'Genetic operators';
gp.operators.mutation.p_mutate = 0.14;    
gp.operators.crossover.p_cross = 0.84;    
gp.operators.directrepro.p_direct = 0.02; 
gp.operators.mutation.mutate_par = [0.9 0.05 0.05 0 0 0];
gp.operators.mutation.gaussian.std_dev = 0.1;  %for mutate_type 3 (constant perturbation): the standard deviation of the Gaussian used.
gp.operators.mutation = orderfields(gp.operators.mutation);
gp.operators = orderfields(gp.operators);

gp.nodes.about = 'Node configuration';
gp.nodes.functions.about = 'Function nodes';
gp.nodes.functions.name = {'times','minus','plus','add3','mult3'};
gp.nodes.functions.arity = [];
gp.nodes.functions.active = [];
gp.nodes.functions = orderfields(gp.nodes.functions);

gp.nodes.const.about = 'Ephemeral random constants';     
gp.nodes.const.num_dec_places = 4;  
gp.nodes.const.range = [-10 10];     %ERC range        
gp.nodes.const.p_ERC = 0.1; %probability of generating an ERC when creating a leaf node
gp.nodes.const.p_int = 0; %probability of generating an integer ERC
gp.nodes.const = orderfields(gp.nodes.const);

gp.nodes.inputs.num_inp = []; 
gp.nodes.inputs.names = cell(0);
gp.nodes.output.name = 'y';
gp.nodes = orderfields(gp.nodes);

%meta data for post-processing
gp.info.filtered = false; %true if population was filtered with GPMODELFILTER
gp.info.lastFilter = []; %the last GPMODELFILTER to be applied
gp.info.merged = 0;  %true if this population is the result of merged independent runs
gp.info.mergedPopSizes = []; %a list of the population sizes that were merged to create the current one
gp.info.duplicatesRemoved = false;
gp.info.version = 'GPTIPS 2 (22.5.15)';

%genes
gp.genes.about = 'Multigene';
gp.genes.multigene = true;                                                                   
gp.genes.max_genes = 4;                                      
gp.genes.operators.p_cross_hi = 0.2;    %probability of high level crossover
gp.genes.operators.hi_cross_rate = 0.5; %probability of any given gene being selected during high level crossover
gp.genes = orderfields(gp.genes);

gp = orderfields(gp);