function res = priv_incrementalSingleBranch(R,analyzer,verbose)
% priv_incrementalSingleBranch - incremental STL verification using a reachable set with a single branch
%
% Syntax:
%    res = priv_incrementalSingleBranch(R,phi,verbose)
%
% Inputs:
%    R - reachable set with just one branch
%    analyzer - onlineReachSetAnalyzer object
%    verbose - whether to print additional information
%
% Outputs:
%    res - four-valued verdict for the current branch
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: priv_modelCheckingIncremental

% Authors:       Florian Lercher
% Written:       15-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that the reachable set has only one branch
if length(R) > 1
    throw(CORAerror('CORA:wrongValue', 'first', 'reachable set with a single branch'))
end

timeInterval = R.timeInterval;
for j = 1:length(timeInterval.time)
    % observe the time interval solution
    t = timeInterval.time{j};
    set = timeInterval.set{j};

    % create STL interval from t
    lb = infimum(t);
    ub = supremum(t);
    % if we are not at the end, excluding the upper time bound is fine,
    % because the time point solution is included in the next interval solution anyway
    rc = j == length(timeInterval.time);
    int = stlInterval(lb,ub,true,rc);
    
    % observe the reachable set
    analyzer.observeSet(set,int,R.loc);
    
    % obtain the current verdict and check if it is conclusive
    v = analyzer.getVerdict();
    if v ~= fourValued.Inconclusive
        res = v;
        if verbose
            disp(['Stop reachability analysis early at time ',num2str(ub)]);
        end
        return;
    end
end

% make sure that all observations are propagated
analyzer.forcePropagation();

res = analyzer.getVerdict();

% ------------------------------ END OF CODE ------------------------------
