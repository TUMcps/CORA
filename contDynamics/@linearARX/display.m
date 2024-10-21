function display(linARX)
% display - Displays a linearARX object on the command window
%
% Syntax:
%    display(linARX)
%
% Inputs:
%    linARX - linearARX object
%
% Outputs:
%    ---
%
% Example:
%    dt = 0.1;
%    A_bar = {[-0.4 0.6; 0.6 -0.4];[0.1 0; 0.2 -0.5]};
%    B_bar = {[0; 0];[0.3; -0.7];[0.1; 0]};
%    linARX = linearARX(A_bar,B_bar,dt)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       02-February-2023 
% Last update:   ---
% Last revision: ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% disp input if necessary
dispInput(inputname(1))

%display parent object
display@contDynamics(linARX);

%display type
disp("Type: Linear discrete-time ARX system");

% display sampling time
disp("Sampling time: " + linARX.dt);

% state equation
disp("y(k) = sum_{i=1}^p A_bar{i} y(k-i) + sum_{i=1}^{p+1} B_bar{i} u(k-i+1)");

% display dimension
disp("Dimension:");
displayMatrixVector(linARX.n_p,"p");
dim_display = min(linARX.n_p, 4);

% display output parameters
disp("Output parameters:");
for i = 1:dim_display
    displayMatrixVector(linARX.A_bar{i,1},sprintf("A_bar%d",i));
end

% display input parameters
disp("Input parameters:");
for i = 1:dim_display+1
    displayMatrixVector(linARX.B_bar{i,1},sprintf("B_bar%d",i));
end

% ------------------------------ END OF CODE ------------------------------
