function res = test_isFuncLinear
% test_isFuncLinear - unit test function for automated check if a function
%    is linear
%
% Syntax:
%    res = test_isFuncLinear()
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       20-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true, wait for failure
res = true;

% for comments:
%    n ... state dimension,
%    m ... input dimension,
%    o ... number of parameters,
%    r ... output dimension.

% 1. One input variable: x

% n = 1, r = 1, linear
f = @(x) x(1);
isLin = isFuncLinear(f);
if ~isLin
    res = false;
end

% n = 3, r = 1, nonlinear
f = @(x) x(1) + x(2)*x(3);
isLin = isFuncLinear(f);
if isLin
    res = false;
end

% n = 3, r = 2, linear + nonlinear
f = @(x) [x(1); x(2)*x(3)];
isLin = isFuncLinear(f);
if ~isLin(1) || isLin(2)
    res = false;
end

% n = 5, r = 2, all linear, matrix multiplication
M = rand(2,5);
f = @(x) M * x(1:5);
isLin = isFuncLinear(f);
if ~all(isLin)
    res = false;
end

% n = 5, r = 2, all linear, matrix multiplication with zeros
M = [0 1 0 0 0; 0 0 1 0 0];
f = @(x) M * x(1:5);
isLin = isFuncLinear(f);
if ~all(isLin)
    res = false;
end


% 2. Two input variables: x, u

% n = 1, m = 1, r = 1, linear
f = @(x,u) x(1) + u(1);
isLin = isFuncLinear(f);
if ~isLin
    res = false;
end

% n = 1, m = 1, r = 1, nonlinear (state)
f = @(x,u) x(1)^2 + u(1);
isLin = isFuncLinear(f);
if isLin
    res = false;
end

% n = 1, m = 1, r = 1, nonlinear (input)
f = @(x,u) x(1) + u(1)^2;
isLin = isFuncLinear(f);
if isLin
    res = false;
end

% n = 4, m = 2, r = 3, linear + nonlinear (state + input)
f = @(x,u) [x(1) + u(2); x(4)^2; exp(u(1))];
isLin = isFuncLinear(f);
if ~isLin(1) || any(isLin(2:3))
    res = false;
end


% 3. Three input variables: x, u, p

% n = 4, m = 2, o = 1, r = 4, linear + nonlinear (state + input + parameters)
f = @(x,u,p) [x(1) + u(2); x(4)^2; p(1)^3; exp(u(1))];
isLin = isFuncLinear(f);
if ~isLin(1) || any(isLin(2:4))
    res = false;
end


% previously failed cases
f = @(in1) reshape( ...
    [1.0, ...
    0.0, ...
    -(in1(4,:).*sin(in1(1,:)).*(9.0./5.0)-in1(3,:).*cos(in1(1,:)).*sin(in1(1,:)).*(8.0./5.0))./(cos(in1(1,:)).^2+1.0)+cos(in1(1,:)).*sin(in1(1,:)).*1.0./(cos(in1(1,:)).^2+1.0).^2.*(in1(4,:).*cos(in1(1,:)).*(9.0./5.0)-in1(3,:).*(cos(in1(1,:)).^2.*(4.0./5.0)-1.0)).*2.0, ...
    -(in1(3,:).*sin(in1(1,:)).*(9.0./5.0)+in1(4,:).*cos(in1(1,:)).*sin(in1(1,:)).*2.0)./(cos(in1(1,:)).^2+1.0)+cos(in1(1,:)).*sin(in1(1,:)).*(in1(4,:).*(cos(in1(1,:)).^2-4.0./5.0)+in1(3,:).*cos(in1(1,:)).*(9.0./5.0)).*1.0./(cos(in1(1,:)).^2+1.0).^2.*2.0, ...
    0.0, ...
    1.0, ...
    0.0, ...
    0.0, ...
    0.0, ...
    0.0, ...
    -(cos(in1(1,:)).^2.*(4.0./5.0)-1.0)./(cos(in1(1,:)).^2+1.0), ...
    (cos(in1(1,:)).*(9.0./5.0))./(cos(in1(1,:)).^2+1.0), ...
    0.0, ...
    0.0, ...
    (cos(in1(1,:)).*(9.0./5.0))./(cos(in1(1,:)).^2+1.0), ...
    (cos(in1(1,:)).^2-4.0./5.0)./(cos(in1(1,:)).^2+1.0)], ...
    [4,4]);
isLin = isFuncLinear(f);
true_isLin = [true,true,true,true;
            true,true,true,true;
            false,true,false,false;
            false,true,false,false];
if ~all(all(isLin == true_isLin))
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
