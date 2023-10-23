function display(trans)
% display - Displays a transition object
%
% Syntax:
%    display(trans)
%
% Inputs:
%    trans - transition object
%
% Outputs:
%    ---
%
% Example:
%    guard = conHyperplane([1 0],0);
%    reset = struct('A',eye(2),'c',zeros(2,1));
%    target = 1;
%    transition(guard,reset,target)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       20-September-2007
% Last update:   18-June-2022 (MW, improved formatting, empty cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf(newline);

disp([inputname(1), ' =']);

fprintf(newline);

% if transition array given...
if length(trans) > 1

    disp("  " + length(trans) + "x1 transition array");

else

    % display guard set
    if (isnumeric(trans.guard) && isempty(trans.guard)) ...
            || (~isnumeric(trans.guard) && ~representsa_(trans.guard,'emptySet',eps))
        disp("Guard set:");
        display(trans.guard);
    else
        disp("Guard set: (none)");
    end
    
    fprintf(newline);
    
    % display reset function
    if isempty(fields(trans.reset))
        % no reset function (empty transition)
        disp("Reset function: (none)");

    elseif isfield(trans.reset,'A')
        % linear reset function
        if trans.reset.hasInput
            disp("Reset function: Ax + Bu + c");
        else
            disp("Reset function: Ax + c");
        end
        
        % display state matrix
        displayMatrixVector(trans.reset.A,"A");
        
        % display input matrix
        if trans.reset.hasInput
            displayMatrixVector(trans.reset.B,"B");
        end
        
        % display constant offset
        displayMatrixVector(trans.reset.c,"c");
        
    elseif isfield(trans.reset,'f')
        % nonlinear reset function
        disp("Reset function: nonlinear");
        
        % use try-block to avoid annoying error messages
        try
            % create symbolic variables
            sys = nonlinearSys(trans.reset.f,trans.reset.stateDim,...
                trans.reset.inputDim);
            vars = symVariables(sys);
            
            % insert symbolic variables into the system equations
            if trans.reset.hasInput
                f = trans.reset.f([vars.x;vars.u]);
            else
                f = trans.reset.f(vars.x);
            end
            % display equations
            for i=1:length(f)
                disp(['  f(',num2str(i),') = ',char(f(i))]);
            end
            fprintf(newline);
        end
    end
    
    % display target
    if isempty(trans.target)
        % empty transition object
        disp("Target location: (none)");
    elseif length(trans.target) == 1
        % ordinary transition (e.g., in hybrid automaton)
        disp("Target location: " + trans.target);
    else
        % parallel hybrid automaton
        disp("Target locations: " + strjoin(string(trans.target),", "));
    end

    fprintf(newline);

    % display synchronization label
    if ~isempty(trans.syncLabel)
        disp("Synchronization label: '" + trans.syncLabel + "'");
    else
        disp("Synchronization label: (none)");
    end

end

fprintf(newline);

% ------------------------------ END OF CODE ------------------------------
