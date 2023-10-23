function vertices = boxPlaneIntersect(hyperbox,hyperplane,beta)
% boxPlaneIntersect - implementation of the intersection algorithms called
%    ExploreBN (Alg.2) and FindSolutions (Alg.3) from [1] calculating the
%    intersection points of an n-dimensional hyperbox with an
%    (n-1)-dimensional hyperplane
%
% Syntax:
%    vertices = boxPlaneIntersect(hyperbox,hyperplane,beta)
%
% Inputs:
%    hyperbox   - (n,2) with lower(:,1)/upper(:,2) bounds
%    hyperplane - val defined by constraining value alpha
%    beta       - scaling vector for hyperplane
%
% Outputs:
%    vertices - vertices of the intersection
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] C. Lara, J. Flores, F. Calderon.
%       "On the Hyperbox - Hyperplane Intersection Problem"

% Authors:       Mark Wetzlinger
% Written:       08-October-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    vertices = [];

%   check if inputs are valid
%   mainly the equation: Aeq * x == b
%   e.g. (1 -0.5)*(x y)' == 2

%   alpha has to be a scalar
    if size(hyperplane) == 1
        alpha = hyperplane;
    else
        throw(CORAerror('CORA:specialError',...
            'The hyperplane has to be defined by a scalar.'));
    end
    
%   beta has to be a line vector (not a matrix)    
    if size(beta,1) > 1
        throw(CORAerror('CORA:specialError','Beta has to be a vector.'));
    end
    
    B = hyperbox;
    offset = 0;
    do_scaling = 0;
    do_translation = 0;

%   1. check for scaling/translation
%   scaling always first!
    dimensionality = size(hyperbox,1);
%   check for scaling (beta must be 1 at all entries for no scaling)
    if ~isequal(beta,ones(1,dimensionality))
        do_scaling = 1;
        B = aux_boxScaling(B,beta);
    end
%   check for translation (lower bound must be zero at all entries)
    if ~isequal(B(:,1),zeros(dimensionality,1))
        do_translation = 1;
        offset = B(:,1)';
        [B,alpha] = aux_boxTranslation(B,alpha);
    end
    
%   2. get vertices (dim = dim of hyperbox)
%   reminder: now ordered hyperbox located at origin
%   hyperplane given as sum(i,n) x_i = alpha
    vertices = aux_exploreBN(B, alpha);

%   3. scale/translate back if needed
%   first translation
    if do_translation == 1
        vertices = aux_translateSolution(vertices,offset,dimensionality);
    end
%   then scaling
    if do_scaling == 1
        vertices = aux_scaleSolution(vertices,beta,dimensionality);
    end
end


% Auxiliary functions -----------------------------------------------------

function L = aux_exploreBN(B, alpha)
%% pseudocode from the paper:
%   INPUT:  ordered hyperbox B located at the origin with the
%           coordinates of the upper corner (b_1,...,b_n)
%           hyperplane defined by alpha
%   OUTPUT: solutions for the intersection of B and alpha
%           adjacent or at the vertices V(s)
%   L <- {}, F <- {(*,...,*)}
%   while F != {} do
%       v <- pop(F)
%       j <- number of defined variables of v
%       w <- (v_1,...,v_j,0,...,0)
%       compute maximum k <= n s.t. sum(i=1,j) vi + sum(i=j+1,k)b_i <= alpha
%       for all i elem {(j+1,...,k)} w_i <- b_i
%       L <- L u aux_findSolutions(w,B,a)
%       if k < n then
%           v' <- (v_1,...,v_j,b_j+1,...,b_k,*,...,*)
%           Z <- P(v,v')
%           for all z elem Z\{v'} push(F,z)
%       end if
%   end while
%   return L

    partZ = 0;

    dim = size(B,1);
    numEdges = dim*2^(dim-1);
    idxL = 1;
% 1 initialization of L and F
    L = zeros(numEdges,dim);
    F = zeros(1,dim) - 1;
    idxF = 1;
% 2 while loop until F is empty
    while (idxF <= size(F,1))
% 3     pop element of F into v
        v = F(idxF,:);
% 4     get number of defined variables of v
        j = find(v == -1,1);
        if (isempty(j))
            j = dim;
        else
            j = j-1;
        end
% 5     initialize all undefined elements with 0 (-> w)
        w = zeros(1,dim);
        for i=1:j
%           assumption: if v_i = 1, in w there should be b_i
            if v(i) == 1
                w(i) = B(i,2);
            end
%           assumption end
        end
% 6     compute max k <= n
        if j ~= dim
            sum_check = 0;
            for i=1:j
                sum_check = sum_check + v(i);
            end
            i = j+1;
            k = j;
            while i <= dim
                sum_check = sum_check + B(i,2);
                if (sum_check > alpha)
                    break;
                end
                k = i;
                i = i+1;
            end
% 7     put b_i at all w_i positions (i from j+1 to k)
            for i=j+1:k
                w(i) = B(i,2);
            end
        else
            k = dim;
        end
% 8     add new solutions to L
%       additional constraint: if sum of w is bigger than alpha, skip
        if sum(w) <= alpha
            newsol = aux_findSolutions(w, B, alpha);
            if ~isempty(newsol)
                numSol = size(newsol,1);
                L(idxL:idxL+numSol-1,:) = newsol;
                idxL = idxL + numSol;
            end
        end
% 9     check for partition
        if partZ == 0 && k < dim
% 10        define new vertex v'
            v_star = zeros(1,dim)-1;
            for i=1:j
                v_star(i) = v(j);
            end
            for i=j+1:k
%               in algorithm v*_(j+1...k)=b_(j+1...k)
                v_star(i) = 1;
            end
% 11        Z ... partition of v and v_star
            Z = aux_createPartition(v, v_star);
% 12        push all elements from Z into F (unless = v')
            F = cat(1,F,Z);
            partZ = 1;
% 13    end if
        end
% 14end while
        idxF = idxF + 1;
    end
% 15return list
% 	remove duplicates
    L = unique(L(1:idxL-1,:), 'rows');
end

function intersectionPoints = aux_findSolutions(v, B, alpha)
%% pseudocode from the paper:
%   INPUT: ordered hyperbox B located at the origin,
%       coordinates of the upper corner (b_1,...,b_n)
%           hyperplane defined by alpha,
%       vertex v = (v_1,...,v_n)
%   OUTPUT: solutions adjacent to or at v
%   L <- {}
%   alpha' <- alpha - sum(i=1,n) v_i
%   if alpha' = 0 then
%       L <- {v}
%   else
%       for all i elem {1,...,n|bi>alpha',vi=0} do
%           t <- v, t_i <- alpha', L <- L u {t}
%       end for
%   end if
%   return L

% 1 initialize intersection point
    intersectionPoints = [];
% 2 calculate new alpha
    alpha_star = alpha;
    for i=1:size(B,1)
        alpha_star = alpha_star - v(i);
    end
% 3 check if v is a solution vertex
    if alpha_star == 0
% 4     v is an intersection point
        intersectionPoints = v;
% 5 v is not an intersection point
    else
% 6     loop from 1 to n including check
        for i=1:size(B,1)
            if B(i,2) >= alpha_star && v(i) == 0
% 7         generate intersection point
                t = v;
                t(i) = alpha_star;
                intersectionPoints = cat(1, intersectionPoints, t);
            end
% 8     end for loop
        end
% 9 end if
    end
% 10return list of intersection points
end

%% AUXILIARY FUNCTIONS I
%   scaling and translation

function scaledBox = aux_boxScaling(B,beta)
%% description:
%   auxiliary function for the scaling of the box
%   this is needed since exploreBN and findSolutions
%   only work with ordered hyperboxes at the origin

%% mathematics
%   following section 6.3 Extensions of the paper
%   hyperbox given      B = [0,b_1]x...x[0,b_n]
%   hyperplane given    sum(i=1,n) beta_i*x_i = alpha
%   hyperbox becomes    B = [0,b_1*beta_1]x...x[0,b_n*beta_n]
%   hyperplane becomes  sum(i=1,n) x_i = alpha

    dimensionBox = size(B,1);
    scaledBox = zeros(dimensionBox,2);
    for i=1:dimensionBox
        if beta(i) < 0
            scaledBox(i,1) = B(i,2)*beta(i);
            scaledBox(i,2) = B(i,1)*beta(i);
        else
            scaledBox(i,:) = B(i,:)*beta(i);
        end
    end
end

function [translatedBox,alpha_star] = aux_boxTranslation(B,alpha)
%% description:
%   auxiliary function for the scaling of the box
%   this is needed since exploreBN and findSolutions
%   only work with ordered hyperboxes at the origin

%% mathematics
%   following section 6.3 Extensions of the paper
%   hyperbox given      [a_1,b_1]x...x[a_n,b_n]
%   hyperplane given    sum(i=1,n) x_i = alpha
%   hyperbox becomes    [0,b_1-a_1]x...x[0,b_n-a_n]
%   hyperplane becomes  sum(i=1,n) x_i_star = alpha - sum(i=1,n) a_i

    dimensionBox = size(B,1);
    translatedBox = zeros(dimensionBox,2);
    alpha_star = alpha;
    for i=1:dimensionBox
        translatedBox(i,2) = B(i,2) - B(i,1);
        alpha_star = alpha_star - B(i,1);
    end
end

function scaledSolution = aux_scaleSolution(solution,vector,dimensionality)
%% description:
%   if the original hyperbox was scaled,
%    this scaled solution has to be re-scaled to the original

%% mathematics:
%   scaled solution:    c* = (c_1*, ..., c_n*)
%   original solution:  c  = (c_1*/beta_1, ..., c_n*/beta_n)

%% code:
    numSolutions = size(solution,1);
    scaledSolution = zeros(numSolutions,dimensionality);
    for i=1:numSolutions
        for j=1:dimensionality
            scaledSolution(i,j) = solution(i,j)/vector(j);
        end
    end
end

function translatedSolution = aux_translateSolution(solution,vector,dimensionality)
%% description:
%   if the original hyperbox was translated,
%    the translated solution has to be re-translated to the original

%% mathematics:
%   translated solution:    c* = (c_1*, ..., c_n*)
%   original solution:      c  = (c_1*+a_1, ..., c_n*+a_n)

%% code:
    numSolutions = size(solution,1);
    translatedSolution = zeros(numSolutions,dimensionality);
    for i=1:numSolutions
        translatedSolution(i,:) = solution(i,:) + vector;
    end
end

%% AUXILIARY FUNCTIONS II
%   manipulating stacks

function stack = aux_fillStack(stack,curr_dim,total_dim,start_index,last_index)
%   recursive function for filling the stack with 0, 1, -1 (for '*')
    if (curr_dim <= total_dim)
        for j=0:(start_index-last_index-1)
            stack(start_index+2*j,:)   = stack(last_index+j,:);
            stack(start_index+2*j,curr_dim)   = 0;
            stack(start_index+2*j+1,:) = stack(last_index+j,:);
            stack(start_index+2*j+1,curr_dim) = 1;
        end
        new_start = start_index+2*j+1+1;
        new_last = start_index;
        % recursive call
        stack = aux_fillStack(stack,curr_dim+1,total_dim,new_start,new_last);
    end
end

function partition = aux_createPartition(v, v_star)
%% description:
%   generates a partition of the two s-faces v and v'
%   v  =    (v_1, ..., v_j, *, ..., ..., *)
%   v' =    (v_1, ..., ..., v_k, *, ..., *)
%   with j < k (note: v_1 to v_j equal in both v and v')

%% mathematics:
%   following Lemma 6.2
%   the partition encloses the set V(v')
%   and all opposites to v' at indices j+1 to k

%% code:
%   initialize partition
    partition = [];
    dimension = size(v_star,2);
%   use fillStack function for the creation of V(v')
    numdef_v_star = find(v_star == -1,1)-1;
    v_star_set = zeros(2^(dimension-numdef_v_star+1)-1,dimension);
%   v_star needed for copying (to be deleted later)
    v_star_set(1,:) = v_star;
    v_star_set = aux_fillStack(v_star_set,numdef_v_star+1,dimension,2,1);
%   beware that V(v') has no more undefined elements
    v_star_set = v_star_set(2^(dimension-numdef_v_star):end,:);
    partition = cat(1,partition,v_star_set);
    
%   get all from i=j+1 to k (between v and v')
    numdef_v = find(v == -1,1)-1;
    new_v = v;
    for i=numdef_v+1:numdef_v_star
        if v_star(i) == 1
            new_v(i) = 0;
            % added 3
            % init v_set with length 2^(dimension-(i+1))-1
            v_set = zeros(2^(dimension-(i-1))-1,dimension);
            v_set(1,:) = new_v;
            v_set = aux_fillStack(v_set,i+1,dimension,2,1);
            partition = cat(1,v_set(2^(dimension-i):end,:),partition);
            new_v(i) = 1;
        elseif v_star(i) == 0
            new_v(i) = 1;
            v_set = zeros(2^(dimension-(i+1))-1,dimension);
            v_set(1,:) = new_v;
            v_set = aux_fillStack(v_set,i+1,dimension,2,1);
            partition = cat(1,v_set(2^(dimension-i):end,:),partition);
            new_v(i) = 0;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
