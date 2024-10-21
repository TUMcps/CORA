function display(SpS)
% display - Displays the properties of a spectraShadow object (center
%    vector, generator matrix, and coefficient matrix) on the command
%    window
%
% Syntax:
%    display(SpS)
%
% Inputs:
%    SpS - spectraShadow object
%
% Outputs:
%    ---
%
% Example: 
%    SpS = spectraShadow([eye(3) eye(3)]);
%    display(SpS);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Adrian Kulmburg
% Written:       01-August-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isemptyobject(SpS)
    
    dispEmptyObj(SpS,inputname(1));
    
else
    
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
    
    %display dimension
    display@contSet(SpS);
    fprintf(newline);
    
    %display center
    disp('c: ');
    disp(full(SpS.c));

    %display generators
    G = full(SpS.G);
    displayGenerators(G,DISPLAYDIM_MAX,'G');
    
    %display coefficient matrix
    A = full(SpS.A);
    displayMatrixVector(A,'A');
    
    disp(" ");
    % display hidden properties
    disp(['Bounded?              ' aux_prop2string(SpS.bounded)]);
    disp(['Empty set?            ' aux_prop2string(SpS.emptySet)]);
    disp(['Full-dimensional set? ' aux_prop2string(SpS.fullDim)]);
    
end

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
