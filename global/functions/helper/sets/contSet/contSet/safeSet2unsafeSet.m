function F = safeSet2unsafeSet(S)
% safeSet2unsafeSet - convert the union of safe sets to an equivalent
%    representation as the union of unsafe sets
%
% Syntax:
%    F = safeSet2unsafeSet(S)
%
% Inputs:
%    S - cell-array storing the safe sets
%
% Outputs:
%    F - cell-array storing the equivalent unsafe sets represented as
%        objects of class polytope
%
% Example: 
%    S{1} = polytope([0 3 3;3 3 0]);
%    S{2} = interval([2;2],[6;5]);
%  
%    F = safeSet2unsafeSet(S);
%  
%    figure; hold on;
%    xlim([-4,10]); ylim([-4,10]);
%    for i = 1:length(F)
%        plot(F{i},[1,2],'FaceColor',CORAcolor("CORA:unsafe"),'EdgeColor','k');
%    end
%    plot(S{1},[1,2],'FaceColor',CORAcolor("CORA:safe"),'EdgeColor','k');
%    plot(S{2},[1,2],'FaceColor',CORAcolor("CORA:safe"),'EdgeColor','k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification, vnnlib2cora

% Authors:       Niklas Kochdumper
% Written:       23-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% represent first safe set by the union of unsafe sets
F = aux_getUnsafeSets(S{1});
nrTotalSets = length(F);

% loop over all safe sets
for i = 2:length(S)
   
    % represent current safe set by the union of unsafe sets
    F_i = aux_getUnsafeSets(S{i});
    nrAddSets = length(F_i);
    
    % compute the intersection with the previous unsafe sets
    F_ = cell(nrTotalSets*nrAddSets,1);
    for j = 1:nrTotalSets
        for k = 1:nrAddSets
            F_{(j-1)*nrAddSets+k,1} = F{j} & F_i{k}; 
        end
    end
    
    % remove empty polytopes
    F = F_(~cellfun(@(P) representsa_(P,'emptySet',eps),F_));
    nrTotalSets = length(F);
end

end


% Auxiliary functions -----------------------------------------------------

function F = aux_getUnsafeSets(S)
% represent the safe set S as a union of unsafe sets 

    % convert to polytope
    P = polytope(S);
    % ensure that halfspace representation is computed
    constraints(P);
    
    % loop over all polytope halfspaces and invert them
    nrCon = length(P.b);
    F = cell(nrCon,1);
    
    for i = 1:nrCon
        F{i} = ~polytope(P.A(i,:),P.b(i)); 
    end
end

% ------------------------------ END OF CODE ------------------------------
