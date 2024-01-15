function display(zB)
% display - Displays the properties of a zonoBundle object (center and 
%    generator matrix of each zonotope) on the command window
%
% Syntax:
%    display(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    ---
%
% Example: 
%    Z1 = zonotope(zeros(2,1),[1 0.5; -0.2 1]);
%    Z2 = zonotope(ones(2,1),[1 -0.5; 0.2 1]);
%    zB = zonoBundle({Z1,Z2});
%    display(zB)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       09-November-2010
% Last update:   02-May-2020 (MW, add empty case)
%                09-June-2020 (MW, remove dependency from zonotope/display)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special cases
if representsa(zB,'emptySet')
    dispEmptySet(zB,inputname(1));
    return
elseif representsa(zB,'fullspace')
    dispRn(zB,inputname(1));
    return
end


fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

% display id, dimension
display@contSet(zB);
fprintf(newline);

% cap number of generators
maxGens = 10;

% display each zonotope
for i=1:zB.parallelSets
    
    disp(['zonotope ',num2str(i),':',newline]);
    
    %display center
    disp('c: ');
    disp(center(zB.Z{i}));

    %display generators
    G = generators(zB.Z{i});
    displayGenerators(G,maxGens,'G');
    
end

% ------------------------------ END OF CODE ------------------------------
