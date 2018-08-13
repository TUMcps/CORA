function res = test_taylm_reexpand(~)
% test_taylm_reexpand - unit_test_function that tests the re-expansion of
%                       a taylor model at the center of a new domain
%
% Syntax:  
%    res = test_taylm_reexpand( ~ )
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


% Test 1: 1D-function -----------------------------------------------------

% create an initial taylor model
int = interval(4,5);
t = taylm(int,10,'x');
c = (rand(4,1)-ones(4,1))*3;
f = @(t) c(1) * t^3 - c(2)*t^5 + c(3)*t - c(4);
T = f(t);

% caclulate the taylor model for the new domain
T_ = reexpand(T,interval(0,2));

% compare the two talyor models
x = -1:0.01:1;
y = zeros(size(x));
y_ = zeros(size(x));

for i = 1:length(x)
    y(i) = evalTaylm(T,x(i));
    y_(i) = evalTaylm(T_,x(i)-1);
end

if any(abs(y-y_) > 1e-10)
    error('Test 1 failed');
end

% % plot the function
% hold on
% plot(x,y,'r');
% plot(x,y_,'b');



% Test 2: 2D-function -----------------------------------------------------

% create an initial taylor model
int_x = interval(4,5);
int_y = interval(1,2);
tx = taylm(int_x,10,'x');
ty = taylm(int_y,10,'y');
c = (rand(4,1)-ones(4,1))*3;
f = @(x,y) c(1) * x^3 * y^2 - c(2)*x^5 * y + c(3)*x * y^3 - c(4);
T = f(tx,ty);

% caclulate the taylor model for the new domain
T_ = reexpand(T,interval([-0.5;-1.5],[1.5;0.5]));

% compare the two talyor models
x1 = -1:0.1:1;
y1 = -1:0.1:1;
[X1,Y1] = meshgrid(x1,y1);

x2 = -1.5:0.1:0.5;
y2 = -0.5:0.1:1.5;
[X2,Y2] = meshgrid(x2,y2);

Z1 = zeros(size(X1));
Z2 = zeros(size(Y1));

for i = 1:size(X1,1)
    for j = 1:size(X1,2)
        x1 = [X1(i,j);Y1(i,j)];
        x2 = [X2(i,j);Y2(i,j)];
        Z1(i,j) = evalTaylm(T,x1);
        Z2(i,j) = evalTaylm(T_,x2);
    end
end

if any(any(abs(Z1-Z2) > 1e-10))
    error('Test 2 failed');
end

% % plot the function
% hold on
% grid on
% 
% s1 = surf(X1,Y1,Z1);
% set(s1,'EdgeColor','none');
% set(s1,'FaceColor','r');
% 
% s2 = surf(X1,Y1,Z2);
% set(s2,'EdgeColor','none');
% set(s2,'FaceColor','b');



% Test 3: 3D-function -----------------------------------------------------

% create an initial taylor model
int_x = interval(4,5);
int_y = interval(1,2);
int_z = interval(3,4);
tx = taylm(int_x,10,'x');
ty = taylm(int_y,10,'y');
tz = taylm(int_z,10,'z');
c = (rand(4,1)-ones(4,1))*3;
f = @(x,y,z) c(1) * x^3 * y^2 * z^4 - c(2)*x^5 * y * z + c(3)*x * y^3 * z^3- c(4);
T = f(tx,ty,tz);

% caclulate the taylor model for the new domain
T_ = reexpand(T,interval([-0.5;-1.5;-1.5],[1.5;0.5;0.5]));

% compare the two talyor models
x1 = -1:0.1:1;
y1 = -1:0.1:1;
z1 = -1:0.1:1;
[X1,Y1,Z1] = meshgrid(x1,y1,z1);

x2 = -1.5:0.1:0.5;
y2 = -0.5:0.1:1.5;
z2 = -0.5:0.1:1.5;
[X2,Y2,Z2] = meshgrid(x2,y2,z2);

diffMax = 0;

for i = 1:size(X1,1)
    for j = 1:size(X1,2)
        for k = 1:size(X1,3)
            x1 = [X1(i,j,k);Y1(i,j,k);Z1(i,j,k)];
            x2 = [X2(i,j,k);Y2(i,j,k);Z2(i,j,k)];
            t1 = evalTaylm(T,x1);
            t2 = evalTaylm(T_,x2);

            diffMax = max(abs(t1-t2),diffMax);
        end
    end
end

if diffMax > 1e-8
    error('Test 3 failed');
end

% If this code is reached, then there was no error and the test is
% considered successful
res = 1;

end



% Auxiliary Functions -----------------------------------------------------

function res = evalTaylm(obj,x)

    mon = obj.monomials;
    coeff = obj.coefficients;
    
    res = 0;
    
    for i = 1:length(coeff)
        e = mon(i,2:end);
        res = res + coeff(i) * prod(x.^(e'));
    end

end



%------------- END OF CODE --------------
