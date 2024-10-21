function display(nlnARX)
% display - displays a nonlinearARX object on the command window
%
% Syntax:
%    display(nlnARX)
%
% Inputs:
%    nlnARX - nonlinearARX object
%
% Outputs:
%    -
%
% Example:
%    f = @(y,u) [y(1,1) + u(1,1) - y(2,1); y(3,1) + u(2,1)*cos(y(1,1)); y(5,1) + u(4,1)*sin(y(1,1))];
%    dt = 0.25;
%    sys = nonlinearARX(f,dt,3,2,2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow
% Written:       04-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% disp input if necessary
dispInput(inputname(1))

%display name and id
disp("Continuous dynamics: '" + nlnARX.name + "'");

% display number of inputs
disp("  dimension of inputs: " + nlnARX.nrOfInputs);

% display number of outputs
disp("  dimension of outputs: " + nlnARX.nrOfOutputs);

% display number of past time steps
disp("  number of past time steps: " + nlnARX.n_p);

fprintf(newline);

%display type
disp('Type: Nonlinear ARX system');

% display sampling time
disp("Sampling time: " + nlnARX.dt);

%create symbolic variables
vars.u = sym('u', 0);
vars.y = sym('y', 0);
for i_t=0:nlnARX.n_p
    if i_t > 0
        for i_dim=1:nlnARX.nrOfOutputs
            command=['vars.y', '(', num2str((i_t-1)*nlnARX.nrOfOutputs+i_dim) ' ,1)=sym(''y',num2str(i_dim),'_t',num2str(i_t),''');'];
            eval(command);
        end
    end
    for i_dim=1:nlnARX.nrOfInputs
        command=['vars.u', '(', num2str(i_t*nlnARX.nrOfInputs+i_dim) ' ,1)=sym(''u',num2str(i_dim),'_t',num2str(i_t),''');'];
        eval(command);
    end
end

%insert symbolic variables into the system equations
f = nlnARX.mFile(vars.y, vars.u);

%display state space equations
disp('Output equations:')
for i_t=1:length(f)
    f_char = strrep(char(f(i_t)), '_t0', '(k)');
    for i_dim=1:nlnARX.n_p
        f_char = strrep(f_char, ['_t',num2str(i_dim)], ['(k-',num2str(i_dim), ')']);
    end
    disp(['  y',num2str(i_t),'(k) = ',f_char]);
end

fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
