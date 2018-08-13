function res = test_globalOptimizer(~)
% test_globalOptimizer - unit_test_function that tests the verified
%                        computation of bounds of a function on some domain
%
% Syntax:  
%    res = test_globalOptimizer( ~ )
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
tol = 0.001;


% Test 1: 1D-Function -----------------------------------------------------

% function f(x) = 1 + x^5 - x^4, x \in [0,1]
f= @(x) 1 + x^5 - x^4;
x = interval(0,1);

% real interval
intReal = interval(0.91808,1);

% calculate the bounds
[int,xMin,domMin,xMax,domMax] = globVerBounds(f,x,tol);

% check if the computed bounds are all over-approximative
checkOverapproximation(intReal,int);

% check if the computed bounds for 'linQuad'-optimization are close enought
% to the real boundaries
checkDiffToOpt(intReal,int,tol);



% Test 2: Beale-Function (2D) ---------------------------------------------

f=@(x) (1.5 - x(1)*(1-x(2))).^2 + (2.25 - x(1)*(1-x(2)^2))^2 + (2.625 - x(1)*(1-x(2)^3))^2;
x = interval([-4.5;-4.5],[4.5;4.5]);

% calculate the bounds 
[int,xMin,domMin,xMax,domMax] = globVerBounds(f,x,tol);

% % plot the function
% [X1,X2] = meshgrid(-4.5:0.1:4.5,-4.5:0.1:4.5);
% Z = zeros(size(X1));
% 
% for i = 1:size(X1,1)
%     for j = 1:size(X1,2)
%         Z(i,j) = f([X1(i,j);X2(i,j)]);
%     end
% end
% 
% hold on
% grid on
% s = surf(X1,X2,Z);
% su = surf(X1,X2,ones(size(Z))*infimum(int));
% set(su,'EdgeColor','none');
% set(su,'FaceColor','r');
% sl = surf(X1,X2,ones(size(Z))*supremum(int));
% set(sl,'EdgeColor','none');
% set(sl,'FaceColor','r');

% real interval
intReal = interval(0,1.818536132812500e+05);

% check if the computed bounds are all over-approximative
checkOverapproximation(intReal,int);

% check if the computed bounds for 'linQuad'-optimization are close enought
% to the real boundaries
checkDiffToOpt(intReal,int,tol);



% Test 3: 2D-function -----------------------------------------------------

f = @(x) -3 + 2*x(1) - x(2) + x(1)^2 + 4 * x(1)*x(2) + 4*x(2)^2;
x = interval([-1;1],[3;4]);

% calculate the bounds
[int,xMin,domMin,xMax,domMax] = globVerBounds(f,x,tol);

% % plot the function
% [X1,X2] = meshgrid(-1:0.1:3,1:0.1:4.);
% Z = zeros(size(X1));
% 
% for i = 1:size(X1,1)
%     for j = 1:size(X1,2)
%         Z(i,j) = f([X1(i,j);X2(i,j)]);
%     end
% end
% 
% hold on
% grid on
% s = surf(X1,X2,Z);
% su = surf(X1,X2,ones(size(Z))*infimum(int));
% set(su,'EdgeColor','none');
% set(su,'FaceColor','r');
% sl = surf(X1,X2,ones(size(Z))*supremum(int));
% set(sl,'EdgeColor','none');
% set(sl,'FaceColor','r');

% real interval
intReal = interval(-5,120);

% check if the computed bounds for 'linQuad'-optimization are close enought
% to the real boundaries
checkDiffToOpt(intReal,int,tol);



% Test 4: Benchmark models ------------------------------------------------

% Define intervals
i_a = interval( -4.5,-0.3); 
i_b = interval(0.4, 0.9);

x = infimum(i_a):0.01:supremum(i_a);
y = infimum(i_b):0.01:supremum(i_b);

% bspline0
z = zeros(size(x));
for i = 1:length(x)
   z(i) = bspline0(x(i)); 
end

int = globVerBounds(@bspline0,i_a,1e-10);

% plot(x,z,'b');
% hold on
% plot(x,infimum(int)*ones(size(x)),'r');
% plot(x,supremum(int)*ones(size(x)),'r');

checkOverapproximationTol(interval(min(z),max(z)),int,1e-10);

% bspline1
z = zeros(size(x));
for i = 1:length(x)
   z(i) = bspline1(x(i)); 
end

int = globVerBounds(@bspline1,i_a,1e-10);

% plot(x,z,'b');
% hold on
% plot(x,infimum(int)*ones(size(x)),'r');
% plot(x,supremum(int)*ones(size(x)),'r');

checkOverapproximationTol(interval(min(z),max(z)),int,1e-10);

% bspline2
z = zeros(size(x));
for i = 1:length(x)
   z(i) = bspline2(x(i)); 
end

int = globVerBounds(@bspline2,i_a,1e-10);

% plot(x,z,'b');
% hold on
% plot(x,infimum(int)*ones(size(x)),'r');
% plot(x,supremum(int)*ones(size(x)),'r');

checkOverapproximationTol(interval(min(z),max(z)),int,1e-10);

% bspline3
z = zeros(size(x));
for i = 1:length(x)
   z(i) = bspline3(x(i)); 
end

int = globVerBounds(@bspline3,i_a,1e-10);

% plot(x,z,'b');
% hold on
% plot(x,infimum(int)*ones(size(x)),'r');
% plot(x,supremum(int)*ones(size(x)),'r');

checkOverapproximationTol(interval(min(z),max(z)),int,1e-10);

% himmilbeau
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));

for i = 1:size(X,1)
    for j = 1:size(X,2)
        Z(i,j) = himmilbeau(X(i,j),Y(i,j));
    end
end

func = @(x) himmilbeau(x(1),x(2));
int = globVerBounds(func,[i_a;i_b],1e-10);

% hold on
% grid on
% s = surf(X,Y,Z);
% 
% s1 = surf(X,Y,infimum(int)*ones(size(X)));
% set(s1,'EdgeColor','none');
% set(s1,'FaceColor','r');
% 
% s2 = surf(X,Y,supremum(int)*ones(size(X)));
% set(s2,'EdgeColor','none');
% set(s2,'FaceColor','r');

checkOverapproximationTol(interval(min(min(Z)),max(max(Z))),int,1e-10);





% Test 5: Lennard-Jones Potential (6D) ------------------------------------

% f = @(x) lennardJonesPotential(x);
% x = interval([0.8;0.4;0.7;0.4;0.2;0.7],[1.2;0.6;1;0.6;0.4;1]);
% 
% [int,xMin,domMin,xMax,domMax] = globVerBounds(f,x,tol,[],'linQuad');





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

function checkOverapproximationTol(intReal,int,tol)
    
    if any(infimum(int) > infimum(intReal)+ones(size(intReal))*tol) || ...
       any(supremum(int) < supremum(intReal)-ones(size(intReal))*tol)
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