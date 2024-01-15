function display(P)
% display - Displays the properties of a polytope object (inequality and
%    equality constraints) on the command window
%
% Syntax:
%    display(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    (to console)
%
% Example: 
%    A = [1 2; -1 2; -2 -2; 1 -2];
%    b = ones(4,1);
%    P = polytope(A,b);
%    display(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Viktor Kotsev, Mark Wetzlinger
% Written:       06-June-2022
% Last update:   01-December-2022 (MW, adapt to other CORA display)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

V = P.V.val;

fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

%display dimension
display@contSet(P);
fprintf(newline);

% display vertex representation
if ~isempty(V)
    disp('Vertex representation:');
    displayMatrixVector(V, 'V');
end

% display inequality constraints
if isempty(P.A)
    disp('Inequality constraints (A*x <= b): (none)');
else
    disp('Inequality constraints (A*x <= b):');
    displayMatrixVector(P.A,'A');
    displayMatrixVector(P.b,'b');
end

% display equality constraints
if isempty(P.Ae)
    disp('Equality constraints (Ae*x = be): (none)');
else
    disp('Equality constraints (Ae*x = be):');
    displayMatrixVector(P.Ae,'Ae');
    displayMatrixVector(P.be,'be');
end

disp(" ");
% display hidden properties
disp(['Bounded?                          ' aux_prop2string(P.bounded)]);
disp(['Empty set?                        ' aux_prop2string(P.emptySet)]);
disp(['Full-dimensional set?             ' aux_prop2string(P.fullDim)]);

disp(['Minimal halfspace representation? ' aux_prop2string(P.minHRep)]);
disp(['Minimal vertex representation?    ' aux_prop2string(P.minVRep)]);

end


% Auxiliary functions -----------------------------------------------------

function res = aux_prop2string(prop)

    if isempty(prop.val)
        res = 'Unknown';
    elseif prop.val
        res = 'true';
    elseif ~prop.val
        res = 'false';
    end

end

% ------------------------------ END OF CODE ------------------------------
