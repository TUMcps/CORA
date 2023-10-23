function res = test_boxPlaneIntersection
% test_boxPlaneIntersection - unit test function for the calculation of
%       vertices of a hyperbox - hyperplane intersection
%
% Syntax:
%    res = test_boxPlaneIntersection
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -
%
% References: 
%   [1] C. Lara, J. Flores, F. Calderon.
%       "On the Hyperbox - Hyperplane Intersection Problem"

% Authors:       Mark Wetzlinger
% Written:       08-October-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
% TEST 1: analytical ------------------------------------------------------

% % type in all values for single check
% dimA = [-5 -4];
% dimB = [ 3  5];
% dim = size(dimA,2);
% hyperbox = [dimA', dimB'];
% beta = [1 2];
% alpha = 2;
% 
% % naive method
% V_naive = boxPlaneIntersectNaive(hyperbox, alpha, beta);
% % proposed method
% V = boxPlaneIntersect(hyperbox, alpha, beta);
% % former method
% A = [diag(ones(dim,1));-diag(ones(dim,1))];
% b = [dimB';-dimA'];
% V_former = lcon2vert(A,b,beta,alpha);
% 
% % compare solutions
% if aux_compareSolutions(V_former, V, V_naive) == 0
%     throw(CORAerror('CORA:testFailed')); 
% end

% TEST 2: random instances ------------------------------------------------

% set minimum and maximum test dimension (for box)
mindim = 3;
maxdim = 5;
% set number of tests per dimension
testperdim = 200;
% no. of tests
testnum = testperdim * (1+maxdim-mindim);
% containers for saving all parameters
allboxesA = zeros(testnum,maxdim);
allboxesB = zeros(testnum,maxdim);
allplanes = zeros(testnum,1);
allbetas = zeros(testnum,maxdim);
% checks = zeros(testnum,1);

% test loop
for d=mindim:maxdim
    for t=1:testperdim
        index = (d-mindim)*testperdim + t;
%       lower and upper bounds of every dimension
        dimA = randi([-5 -3],1,d);
        dimB = randi([3 5],1,d);
        allboxesA(index,1:d) = dimA;
        allboxesB(index,1:d) = dimB;
%       hyperbox(1,1) = a_1, hyperbox(1,2) = b_1
%       hyperbox(2,1) = a_2, ...
        hyperbox = [dimA', dimB'];
%       constraint for hyperplane: sum(i,n) x_i = alpha
        alpha = randi([1 3],1);
        allplanes(index) = alpha;
        beta = randi([1 2],1,d);
        allbetas(index,1:d) = beta;
%       former method
        A = [diag(ones(d,1));-diag(ones(d,1))];
        b = [dimB';-dimA'];
%         V_former = lcon2vert(A,b,beta,alpha);
%       naive method
%         V_naive = boxPlaneIntersectNaive(hyperbox, alpha, beta);
%       proposed method
        V = boxPlaneIntersect(hyperbox, alpha, beta);
%       compare solutions
%         same = aux_compareSolutions(V_former, V_naive, V);
%         checks(index) = same;
    end
end

% misscount = size(checks,1) - nnz(checks);
% if misscount > 0
%     throw(CORAerror('CORA:testFailed'));
% end

res = true;

end


% Auxiliary functions -----------------------------------------------------

function eqty = aux_compareSolutions(sol_1, sol_2, sol_3, TOL)
%% description:
%   checks if all solutions are the same
%   for 2 or 3 lists of solution points
%   TOL optional
%   returns 1 same, 0 if not

%% code:
    eqty = 1;
    if nargin < 4
        TOL = 1e-6;
    end
    if nargin == 3
%   must have same number of intersection points
        if size(sol_2,1) ~= size(sol_3,1)
            eqty = 0;
            return;
        elseif size(sol_1,1) ~= size(sol_2,1)
            eqty = 0;
            return;
        end
%       intersection points must match
        i = 1;
        while i <= size(sol_1,1)
%       always take first one in sol_2 as to-be-checked point
            for j=1:size(sol_2,1)
%           search if same point in sol_3
                if abs(sol_1(i,:) - sol_2(j,:)) < TOL
                    for k=1:size(sol_3,1)
%                   search if same point in former
                        if abs(sol_3(k,:) - sol_2(j,:)) < TOL
                            sol_1 = cat(1,sol_1(1:i-1,:),sol_1(i+1:end,:));
                            sol_2 = cat(1,sol_2(1:j-1,:),sol_2(j+1:end,:));
                            sol_3 = cat(1,sol_3(1:k-1,:),sol_3(k+1:end,:));
                            j = -1;
                            k = -1;
                            break;
                        end
                    end
                    if k == size(sol_3,1)
%                   no matching solution in sol_3 found -> return failure
                        eqty = 0;
                        return;
                    end
                    break;
                end
            end
            if j == size(sol_2,1)
%           no matching solution in sol_2 found -> return failure
                eqty = 0;
                return;
            end
        end
        
%       check if all lists now empty (meaning all had same points)
        if size(sol_1,1) > 0 || size(sol_2,1) > 0 || size(sol_3,1) > 0
            eqty = 0;
        end
    elseif nargin < 3
%   only two lists to be checked
        if size(sol_1,1) ~= size(sol_2,1)
            eqty = 0;
            return;
        end
%       intersection points must match
        i = 1;
        while i <= size(sol_1,1)
%       always take first one in sol_1 as to-be-checked point
            for j=1:size(sol_2,1)
%           search if same point in sol_2
                if abs(sol_1(i,:) - sol_2(j,:)) < TOL
                    sol_1 = cat(1,sol_1(1:i-1,:),sol_1(i+1:end,:));
                    sol_2 = cat(1,sol_2(1:j-1,:),sol_2(j+1:end,:));
                    j = -1;
                    break;
                end
            end
            if j == size(sol_2,1)
%           no matching solution in sol_2 found -> return failure
                eqty = 0;
                return;
            end
        end
        
%       check if all lists now empty (meaning all had same points)
        if size(sol_1,1) > 0 || size(sol_2,1) > 0
            eqty = 0;
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
