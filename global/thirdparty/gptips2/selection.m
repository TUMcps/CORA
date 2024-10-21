function ID = selection(gp)
%SELECTION Selects an individual from the current population.
%
%   Selection is based on fitness and/or complexity using either regular or
%   Pareto tournaments.
%
%   ID = SELECTION(GP) returns the numeric identifier of the selected
%   individual.
%
%   Remarks:
%
%   Selection method is either:
%
%   A) Standard fitness based tournament (with the option of lexicographic
%   selection pressure).
%
%   This is always used when the field GP.SELECTION.TOURNAMENT.P_PARETO
%   equals 0.
%
%   OR
%
%   B) Pareto tournament
%
%   This is used with a probability equal to
%   GP.SELECTION.TOURNAMENT.P_PARETO
%
%   Remarks:
%
%   Method A (standard tournament)
%   ------------------------------
%
%   For a tournament of size N
%
%   1) N individuals are randomly selected from the population with
%   reselection allowed.
%
%   2) The population index of the best individual in the tournament is
%   returned.
%
%   3) If one or more individuals in the tournament have the same best
%   fitness then one of these is selected randomly, unless
%   GP.SELECTION.TOURNAMENT.LEX_PRESSURE is set to true. In this case the
%   one with the lowest complexity is selected (if there is more than one
%   individual with the best fitness AND the same complexity then one of
%   these individuals is randomly selected.)
%
%   OR
%
%   Method B (Pareto tournament)
%   ----------------------------
%
%   This is used probabalistically when GP.SELECTION.TOURNAMENT.P_PARETO > 0
%
%   For a tournament of size N
%
%   1) N individuals are randomly selected from the population with
%   reselection allowed.
%
%   2) These N individuals are then Pareto sorted (sorted using fitness and
%   complexity as competing objectives) to obtain a (rank 1) pareto set of
%   individuals that are non-dominated in terms of both fitness and
%   complexity (this set can be thought of as being the pareto front of the
%   tournament). The pareto sort is based on on the sorting method
%   described on page 184 of:
%
%   "A fast and elitist multiobjective genetic algorithm: NSGA-II" by
%   Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal, T. Meyarivan IEEE
%   Transactions on Evolutionary Computation VOL. 6, NO. 2, APRIL 2002, pp.
%   182-197.
%
%   3) A random individual is then selected from the pareto set and its
%   index IND in the population is returned.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also POPBUILD, INITBUILD

%Method A
if rand >= gp.selection.tournament.p_pareto
    
    %pick N individuals at random
    tour_ind = ceil(rand(gp.selection.tournament.size,1)*gp.runcontrol.pop_size);
    
    %retrieve fitness values of tournament members
    tour_fitness = gp.fitness.values(tour_ind);
    
    %check if max or min problem
    if gp.fitness.minimisation
        bestfitness = min(tour_fitness);
    else
        bestfitness = max(tour_fitness);
    end
    
    %locate indices of best individuals in tournament
    bestfitness_tour_ind = find(tour_fitness == bestfitness);
    number_of_best = numel(bestfitness_tour_ind);
    
    %use plain lexicographic parsimony pressure a la Sean Luke
    %and Livui Panait (Lexicographic Parsimony Pressure, GECCO 2002, pp. 829-836).
    if number_of_best > 1
        
        if gp.selection.tournament.lex_pressure
            
            tcomps = gp.fitness.complexity(tour_ind(bestfitness_tour_ind));
            
            [~,min_ind] = min(tcomps);
            bestfitness_tour_ind = bestfitness_tour_ind(min_ind);
            
        else %otherwise randomly pick one
            
            bestfitness_tour_ind = bestfitness_tour_ind(ceil(rand*number_of_best));
            
        end
    end
    
    ID = tour_ind(bestfitness_tour_ind);
     
else %Method B: Pareto tournament
    
    %pick N individuals at random
    tour_ind = ceil(rand(gp.selection.tournament.size,1)*gp.runcontrol.pop_size);
    
    %retrieve fitness values of tournament members
    tour_fitness = gp.fitness.values(tour_ind);
    
    %retrieve complexities of tournament members
    tour_comp  = gp.fitness.complexity(tour_ind);
    
    %perform fast pareto sort on tournament members
    %and select winner randomly from rank 1 solutions
    if gp.fitness.minimisation
        mo = [tour_fitness tour_comp];
    else
        mo = [ (-1 *tour_fitness) tour_comp];
    end
    
    %remove any 'infs' from consideration and recall if necessary
    infs= isinf(mo(:,1));
    mo(infs,:)=[];
    if isempty(mo)
        ID = selection(gp);
        return;
    end
    
    tour_ind(infs,:)=[];
    
    rank = ndfsort_rank1(mo);
    rank1_tour_ind = tour_ind(rank == 1);
    num_rank1 = numel(rank1_tour_ind);
    ID = rank1_tour_ind(ceil(rand * num_rank1));
end