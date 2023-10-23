function F = safeSet2unsafeSet(S)
% safeSet2unsafeSet - convert the union of safe sets to an equivalent
%                     representation as the union of unsafe sets
%
% Syntax:
%    F = safeSet2unsafeSet(S)
%
% Inputs:
%    S - cell-array storing the safe sets
%
% Outputs:
%    F - cell-array sotring the equivalent unsafe sets represented as
%        objects of class polytope
%
% Example: 
%    S{1} = polytope([0 3 3;3 3 0]);
%    S{2} = interval([2;2],[6;5]);
%  
%    F = safeSet2unsafeSet(S);
% 
%    csafe = CORAcolor("CORA:safe");
%    cunsafe = CORAcolor("CORA:unsafe");
%  
%    figure; hold on;
%    xlim([-4,10]); ylim([-4,10]);
%    for i = 1:length(F)
%       plot(F{i},[1,2],'FaceColor',cunsafe,'EdgeColor','k');
%    end
%    plot(S{1},[1,2],'FaceColor',csafe,'EdgeColor','k');
%    plot(S{2},[1,2],'FaceColor',csafe,'EdgeColor','k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Authors:       Niklas Kochdumper
% Written:       23-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    F = aux_getUnsafeSets(S{1});

    % loop over all safe sets
    for i = 2:length(S)
       
        % represent current safe set by the union of unsafe sets
        Ftemp = aux_getUnsafeSets(S{i});
        
        % compute the intersection with the previous unsafe sets
        F_ = cell(length(F)*length(Ftemp),1);
        
        for j = 1:length(F)
           for k = 1:length(Ftemp)
              F_{(j-1)*length(Ftemp)+k,1} = F{j} & Ftemp{k}; 
           end
        end
        
        % remove empty polytopes
        F = F_(~cellfun(@(P) representsa_(P,'emptySet',eps),F_));
    end
end


% Auxiliary functions -----------------------------------------------------

function F = aux_getUnsafeSets(S)
% represent the safe set S as a union of unsafe sets 

    % convert to polytope
    S = polytope(S);
    
    % loop over all polytope halfspaces and invert them
    F = cell(size(S.A,1),1);
    
    for i = 1:length(F)
       F{i} = polytope(-S.A(i,:),-S.b(i)); 
    end
end

% ------------------------------ END OF CODE ------------------------------
