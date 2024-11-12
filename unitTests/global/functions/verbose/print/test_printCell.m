function res = test_printCell
% test_printCell - unit test function for printCell
%
% Syntax:
%    res = test_printCell()
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

% Authors:       Maximilian Perschl
% Written:       05-November-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty
printCell({})

% scalar
printCell({1})
printCell({"text"})
printCell({1,"text"})
printCell({1,"text",[1;2;3]})
printCell({1,"text",[1 2 3;4 5 6],interval(2,3)})

C = {1,"text",[1 2 3;4 5 6],interval(2,3)};

% row vector
printCell([C,C])
printCell([C,C,C,C,C])

% column vector
printCell([C;C])
printCell([C;C;C;C;C])

% matrix
printCell([C C;C C])

% parameters
printCell(C);
printCell(C,'%4.3f%s');
printCell(C,'high');
printCell(C,'high',true);
printCell(C,'high',false);
fprintf('\n')

% example
C = struct('a',[1 2 3],'b','text');
printStruct(C)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
