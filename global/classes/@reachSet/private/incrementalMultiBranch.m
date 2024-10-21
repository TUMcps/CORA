function [res, contagious] = incrementalMultiBranch(R,analyzer,i,tFinal,verbose)
% incrementalMultiBranch - incremental STL verification using a reachable set with multiple branches
%
% Syntax:
%   [res, contagious] = incrementalMultiBranch(R,phi,i,tFinal,verbose)
%
% Inputs:
%    R - reachable set
%    analyzer - onlineReachSetAnalyzer object
%    i - index of the current branch
%    tFinal - final time step of R
%    verbose - whether to print additional information
%
% Outputs:
%    res - four-valued verdict for the current branch
%    contagious - true if inconclusive results should be propagated to other branches
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: modelCheckingIncremental

% Authors:       Florian Lercher
% Written:       15-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

cs = children(R,i);
if isempty(cs)
    tStartChildren = []; 
else
    [tStartChildren,idx] = sort(arrayfun(@(c) infimum(R(c).timeInterval.time{1}),cs));
    cs = cs(idx);
end

[res,tStop,checkpoints] = aux_single_branch(R(i),analyzer,tStartChildren);
tFinal_branch = supremum(R(i).timeInterval.time{end});
if verbose && res ~= fourValued.Inconclusive && tStop < tFinal_branch
    disp(['Stop reachability branch ',num2str(i),' early at time ',num2str(tStop)]);
end
% if we reached the final time step, an inconclusive verdict is contagious
contagious = tStop == tFinal;

for c = 1:length(cs)
    new_d = checkpoints{c};
    if isempty(new_d)
        new_d = analyzer.copy();
    end
    tStart = infimum(R(cs(c)).timeInterval.time{1});
    % we only need to look at the new branch if it started before we obtained a verdict
    if tStart <= tStop
        [branch_res,branch_cont] = incrementalMultiBranch(R,new_d,cs(c),tFinal,verbose);
        res = aux_combine_results(res,branch_res,contagious,branch_cont);
        contagious = contagious || branch_cont;
    elseif verbose
        disp(['Skip reachability branch ',num2str(cs(c)),' and its children']);
    end
end
end


% Auxiliary functions -----------------------------------------------------

function [res,tStop,checkpoints] = aux_single_branch(R,analyzer,tStartChildren)
    assert(length(R) == 1);
    timeInterval = R.timeInterval;
    checkpoints = cell(size(tStartChildren));
    numChildren = length(tStartChildren);
    nextCheckpoint = 1;

    for j = 1:length(timeInterval.time)
        % observe the time interval solution
        int = timeInterval.time{j};
        lb = infimum(int);
        ub = supremum(int);
        
        makeValid = stlInterval(0,lb,true,false);
        if nextCheckpoint <= numChildren && tStartChildren(nextCheckpoint) <= ub
            % create checkpoint for child branches that start now
            % before adding observations of this branch
            int = stlInterval(lb,tStartChildren(nextCheckpoint),true,false);
            lb = tStartChildren(nextCheckpoint);
            analyzer.observeSet(timeInterval.set{j},int,R.loc,makeValid);

            while nextCheckpoint <= numChildren && tStartChildren(nextCheckpoint) <= ub
                % create checkpoint for child branches that start now
                checkpoints{nextCheckpoint} = analyzer.copy();
                nextCheckpoint = nextCheckpoint + 1;
            end
        end

        % when we are at the final time step, we know that no more observations will be added
        if j == length(timeInterval.time)
            makeValid = stlInterval(0,ub);
            rc = true;
        else
            rc = false;
        end
        analyzer.observeSet(timeInterval.set{j},stlInterval(lb,ub,true,rc),R.loc,makeValid);
        v = analyzer.getVerdict();
        if v ~= fourValued.Inconclusive
            res = v;
            tStop = ub;
            return;
        end
    end

    analyzer.forcePropagation();
    res = analyzer.getVerdict();
    tStop = supremum(timeInterval.time{end});
end

function res = aux_combine_results(base_res,new_res,base_contagious,new_contagious)
    if base_contagious && base_res == fourValued.Inconclusive
        res = base_res;
    elseif new_contagious && new_res == fourValued.Inconclusive
        res = new_res;
    elseif base_res ~= fourValued.Inconclusive && new_res ~= fourValued.Inconclusive
        if base_res == new_res
            res = base_res;
        else
            res = kleene.Unknown;
        end
    elseif base_res ~= fourValued.Inconclusive
        res = base_res;
    else
        res = new_res;
    end
end

% ------------------------------ END OF CODE ------------------------------
