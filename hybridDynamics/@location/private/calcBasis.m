function D = calcBasis(obj,R,Pguard,options)
% calcBasis - calculate suitable orthogonal transformation matrices for the
%             tight overapproximation of the sets R with interval
%
% Syntax:  
%    D = calcBasis(obj,R,Pguard,options)
%
% Inputs:
%    obj - location object
%    R - cell array of reachable sets
%    Pguard - guard set that was intersected by the reachable sets
%    options - struct containing algorithm options
%
% Outputs:
%    D - cell array containing the different transformation matrices
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      22-May-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------



    % initialization
    D = cell(length(options.enclosureEnables),1);
    
    % loop over all selected options
    for i = 1:length(options.enclosureEnables)
        
        switch options.enclosureEnables(i)
           
            % axis aligned bounding box
            case 1
                D{i} = eye(obj.contDynamics.dim);
                
                
                
            % Principal Component Analysis (generators)
            case 2
                
                % extract hyperplane normal vector
                if isa(Pguard,'halfspace')
                   n = get(Pguard,'c'); 
                elseif isa(Pguard,'constrainedHyperplane')
                   n = get(Pguard.h,'c');
                end
                
                % normalize normal vector
                n = n./norm(n);
                
                % concatenate all generators
                G = extractGenerators(R);
                
                % project the generators onto the hyperplane
                if isa(Pguard,'halfspace') || isa(Pguard,'constrainedHyperplane')
                   G = G - n * n'*G;
                end
                
                % calcualte an orthogonal transformation matrix using PCA
                [D{i},~,~] = svd(G); 
                
                
                
            % Principal Component Analysis (center)
            case 3
                
                % extract hyperplane normal vector
                if isa(Pguard,'halfspace')
                   n = get(Pguard,'c'); 
                elseif isa(Pguard,'constrainedHyperplane')
                   n = get(Pguard.h,'c');
                end
                
                % normalize normal vector
                n = n./norm(n);
                
                % concatenate all zonotope mid points
                c = extractMidPoints(R);
                
                % project the generators onto the hyperplane
                if isa(Pguard,'halfspace') || isa(Pguard,'constrainedHyperplane')
                   c = c - n * n'*c;
                end
                
                % calcualte an orthogonal transformation matrix using PCA
                [D{i},~,~] = svd(c); 
            
            
        end        
    end
end



% Auxiliary Functions -----------------------------------------------------

function G = extractGenerators(R)
% extract all generator vectors from the sets that are over-approximated

    G = [];

    % loop over all sets
    for j = 1:length(R)
       if iscell(R{j})              % sets are split
          for k = 1:length(R{j})
              if isa(R{j}{k},'zonotopeBundle')
                    for i = 1:R{j}{k}.parallelSets
                        G = [G,R{j}{k}.Z{i}.Z(:,2:end)];
                    end
              else
                    G = [G,R{j}{k}.Z(:,2:end)];
              end
          end
       else                         % sets are not split
          if isa(R{j},'zonotopeBundle')
                for i = 1:R{j}.parallelSets
                    G = [G,R{j}.Z{i}.Z(:,2:end)];
                end
          else
                G = [G,R{j}.Z(:,2:end)];
          end
       end
    end
end


function c = extractMidPoints(R)
% extract all generator vectors from the sets that are over-approximated

    c = [];

    % loop over all sets
    for j = 1:length(R)
       if iscell(R{j})              % sets are split
          for k = 1:length(R{j})
              if isa(R{j}{k},'zonotopeBundle')
                    for i = 1:R{j}{k}.parallelSets
                        c = [c,R{j}{k}.Z{i}.Z(:,1)];
                    end
              else
                    c = [c,R{j}{k}.Z(:,1)];
              end
          end
       else                         % sets are not split
          if isa(R{j},'zonotopeBundle')
                for i = 1:R{j}.parallelSets
                    c = [c,R{j}.Z{i}.Z(:,1)];
                end
          else
                c = [c,R{j}.Z(:,1)];
          end
       end
    end
end

%------------- END OF CODE --------------