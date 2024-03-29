function matP = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of two matrix
%    polytopes or a matrix polytope with a matrix
%
% Syntax:
%    matP = plus(summand1,summand2)
%
% Inputs:
%    summand1 - matPolytope object or numerical matrix
%    summand2 - matPolytope object or numerical matrix
%
% Outputs:
%    matP - matrix polytope after Minkowski addition
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       21-June-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%Find a matrix polytope object
[matP,summand] = findClassArg(summand1,summand2,'matPolytope');

%Is summand a matrix polytope?
if isa(summand,'matPolytope')
    %initialize matrix vertices
    for i=1:summand1.verts
        V1(i,:)=mat2vec(summand1.vertex{i});
    end
    for i=1:summand2.verts
        V2(i,:)=mat2vec(summand2.vertex{i});
    end
    %initialize potential vertices
    Vpot=[];
    %Calculate posiible new vertices by adding all combinations
    for j=1:length(V1(:,1))
        for i=1:length(V2(:,1))
            Vpot(end+1,:)=V1(j,:)+V2(i,:);
        end
    end
    %compute convex hull
    try
        opt{1}='QJ';
        K=convhulln(Vpot, opt);

        %rewrite result as a matrix polytope
        matV=[];
        for i=1:length(Vpot(:,1))
            if find(K==i)
                matV{end+1}=vec2mat(Vpot(i,:));
            end
        end
    catch
        %rewrite result as a matrix polytope
        matV=[];
        for i=1:length(Vpot(:,1))
            matV{end+1}=vec2mat(Vpot(i,:));
        end
        disp('convex hull could not be computed')
    end
    
    matP=matPolytope(matV);
    
%is summand a matrix?
elseif isnumeric(summand)
    %Calculate minkowski sum
    for i=1:matP.vertices
        matP.vertex{i} = matP.vertex{i} + summand;
    end
end

% ------------------------------ END OF CODE ------------------------------
