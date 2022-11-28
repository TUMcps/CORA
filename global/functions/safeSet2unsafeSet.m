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
%        objects of class mptPolytope
%
% Example: 
%    S{1} = mptPolytope([0 3 3;3 3 0]');
%    S{2} = interval([2;2],[6;5]);
%
%    F = safeSet2unsafeSet(S);
%
%    figure; hold on;
%    xlim([-4,10]); ylim([-4,10]);
%    for i = 1:length(F)
%       plot(F{i},[1,2],'FaceColor','r','EdgeColor','k');
%    end
%    plot(S{1},[1,2],'FaceColor','g','EdgeColor','k');
%    plot(S{2},[1,2],'FaceColor','g','EdgeColor','k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Author:       Niklas Kochdumper
% Written:      23-November-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    F = getUnsafeSets(S{1});

    % loop over all safe sets
    for i = 2:length(S)
       
        % represent current safe set by the union of unsafe sets
        Ftemp = getUnsafeSets(S{i});
        
        % compute the intersection with the previous unsafe sets
        F_ = cell(length(F)*length(Ftemp),1);
        
        for j = 1:length(F)
           for k = 1:length(Ftemp)
              F_{(j-1)*length(Ftemp)+k,1} = F{j} & Ftemp{k}; 
           end
        end
        
        % remove empty polytopes
        F = F_(~cellfun(@(x) isempty(x),F_));
    end
end


% Auxiliary Functions -----------------------------------------------------

function F = getUnsafeSets(S)
% represent the safe set S as a union of unsafe sets 

    % convert to polytope
    S = mptPolytope(S);
    
    % loop over all polytope halfspaces and invert them
    F = cell(size(S.P.A,1),1);
    
    for i = 1:length(F)
       F{i} = mptPolytope(-S.P.A(i,:),-S.P.b(i)); 
    end
end

%------------- END OF CODE --------------