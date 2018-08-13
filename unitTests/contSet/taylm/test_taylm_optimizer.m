function res = test_taylm_optimizer(~)
% test_taylm_optimizer - unit_test_function that tests all optimization 
%                        techniques that can be used to determine the
%                        bounds of a taylor model on some domain
%
% Syntax:  
%    res = test_taylm_optimizer( ~ )
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Niklas Kochdumper
% Written:      14-April-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


res = 0;


% Test 1: 1D-Function -----------------------------------------------------

% function f(x) = 1 + x^5 - x^4, x \in [0,1]
f= @(x) 1 + x^5 - x^4;
x = interval(0,1);
t = taylm(x,10,'x');

% real interval
intReal = interval(0.91808,1);

% calculate the bounds with different optimizers
t = set(t,'opt_method','int');
T = f(t);
intInt = interval(T);

t = set(t,'opt_method','bnb');
T = f(t);
intBnb = interval(T);

t = set(t,'opt_method','bnbAdv');
T = f(t);
intBnbAdv = interval(T);

t = set(t,'opt_method','linQuad');
T = f(t);
intLinQuad = interval(T);

% check if the computed bounds are all over-approximative
checkOverapproximation(intReal,intInt);
checkOverapproximation(intReal,intBnb);
checkOverapproximation(intReal,intBnbAdv);
checkOverapproximation(intReal,intLinQuad);

% check if the computed bounds for 'linQuad'-optimization are close enought
% to the real boundaries
checkDiffToOpt(intReal,intLinQuad,t.eps);





% Test 2: Beale-Function (2D) ---------------------------------------------

f=@(x1,x2) (1.5 - x1*(1-x2)).^2 + (2.25 - x1*(1-x2^2))^2 + (2.625 - x1*(1-x2^3))^2;

x1 = interval(-4.5,4.5);
x2 = interval(-4.5,4.5);

t1 = taylm(x1,10,'x1');
t2 = taylm(x2,10,'x2');

% plot the function
% [X1,X2] = meshgrid(-4.5:0.1:4.5,-4.5:0.1:4.5);
% for i = 1:size(X1,1)
%     for j = 1:size(X1,2)
%         Z(i,j) = f(X1(i,j),X2(i,j));
%     end
% end
% s = surf(X1,X2,Z);

% real interval
intReal = interval(0,1.818536132812500e+05);

% calculate the bounds with different optimizers
t1 = set(t1,'opt_method','int');
t2 = set(t2,'opt_method','int');
T = f(t1,t2);
intInt = interval(T);

t1 = set(t1,'opt_method','bnb');
t2 = set(t2,'opt_method','bnb');
T = f(t1,t2);
intBnb = interval(T);

t1 = set(t1,'opt_method','bnbAdv');
t2 = set(t2,'opt_method','bnbAdv');
T = f(t1,t2);
intBnbAdv = interval(T);

t1 = set(t1,'opt_method','linQuad');
t2 = set(t2,'opt_method','linQuad');
T = f(t1,t2);
intLinQuad = interval(T);

% check if the computed bounds are all over-approximative
checkOverapproximation(intReal,intInt);
checkOverapproximation(intReal,intBnb);
checkOverapproximation(intReal,intBnbAdv);
checkOverapproximation(intReal,intLinQuad);

% check if the computed bounds for 'linQuad'-optimization are close enought
% to the real boundaries
checkDiffToOpt(intReal,intLinQuad,t1.eps);





% Test 3: 2D-function -----------------------------------------------------

f = @(x1,x2) -3 + 2*x1 - x2 + x1^2 + 4 * x1*x2 + 4*x2^2;

x1 = interval(-1,3);
x2 = interval(1,4);

t1 = taylm(x1,10,'x1');
t2 = taylm(x2,10,'x2');

% plot the function
% [X1,X2] = meshgrid(-1:0.1:3,1:0.1:4.);
% Z = zeros(size(X1));
% 
% for i = 1:size(X1,1)
%     for j = 1:size(X1,2)
%         Z(i,j) = f(X1(i,j),X2(i,j));
%     end
% end
% s = surf(X1,X2,Z);

% real interval
intReal = interval(-5,120);

% calculate the bounds with different optimizers
t1 = set(t1,'opt_method','int');
t2 = set(t2,'opt_method','int');
T = f(t1,t2);
intInt = interval(T);

t1 = set(t1,'opt_method','bnb');
t2 = set(t2,'opt_method','bnb');
T = f(t1,t2);
intBnb = interval(T);

t1 = set(t1,'opt_method','bnbAdv');
t2 = set(t2,'opt_method','bnbAdv');
T = f(t1,t2);
intBnbAdv = interval(T);

t1 = set(t1,'opt_method','linQuad');
t2 = set(t2,'opt_method','linQuad');
T = f(t1,t2);
intLinQuad = interval(T);

% check if the computed bounds are all over-approximative
checkOverapproximation(intReal,intInt);
checkOverapproximation(intReal,intBnb);
checkOverapproximation(intReal,intBnbAdv);
checkOverapproximation(intReal,intLinQuad);

% check if the computed bounds for 'linQuad'-optimization are close enought
% to the real boundaries
checkDiffToOpt(intReal,intLinQuad,t1.eps);





% Test 4: Lennard-Jones Potential (6D) ------------------------------------

% x = interval([0.8;0.4;0.7;0.4;0.2;0.7],[1.2;0.6;1;0.6;0.4;1]);
% tx = taylm(x,10,'x','linQuad');
% 
% T = lennardJonesPotential(tx);
% 
% int = globalOptimizer(T)
% 
% 



res = 1;


end



% Auxiliary Functions -----------------------------------------------------


function f = lennardJonesPotential(x)

    a{1} = {0;0;0};
    a{2} = {x(1);0;0};
    a{3} = {x(2);x(3);0};
    a{4} = {x(4);x(5);x(6)};
    
    V = @(r) 1/r.^12 - 2/r.^6;
    
    f = 6;
    for j = 1:4
        for i = 1:j
            r = sqrt((a{i}{1}-a{j}{1})^2 + (a{i}{2}-a{j}{2})^2 + (a{i}{3}-a{j}{3})^2 );
            f = f + V(r);
        end
    end
end

function checkOverapproximation(intReal,int)
    
    if any(infimum(int) > infimum(intReal)) || any(supremum(int) < supremum(intReal))
        error('test_taylm_optimizer failed!');
    end
end

function checkDiffToOpt(intReal,int,tol)
% check if the computed boundaries are close enough to the real boundaries
% of the function

    if any(abs(infimum(intReal) - infimum(int)) > tol) || any(abs(supremum(int) - supremum(intReal)) > tol)
        error('test_taylm_optimizer failed!');
    end
end

%------------- END OF CODE --------------