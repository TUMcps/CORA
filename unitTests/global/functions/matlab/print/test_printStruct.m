function res = test_printStruct
% test_printStruct - unit test function for printStruct
%
% Syntax:
%    res = test_printStruct()
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

% Authors:       Tobias Ladner
% Written:       31-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% empty
printStruct(struct)
resvec(end+1) = true;

% scalar
printStruct(struct('a',1))
printStruct(struct('b','text'))
printStruct(struct('a',1,'b','text'))
printStruct(struct('a',1,'b','text','c',[1 2 3; 4 5 6]))
resvec(end+1) = true;

S = struct('a',1,'b','text','c',[1 2 3; 4 5 6]);

% row vector
printStruct([S,S])
printStruct([S,S,S,S,S])
resvec(end+1) = true;

% column vector
printStruct([S;S])
printStruct([S;S;S;S;S])
resvec(end+1) = true;

% matrix
printStruct([S S;S S])

% parameters
printStruct(S);
printStruct(S,'%4.3f%s');
printStruct(S,'high');
printStruct(S,'high',true);
printStruct(S,'high',false);
fprintf('\n')
resvec(end+1) = true;

% example
S = struct('a',[1 2 3],'b','text');
printStruct(S)
resvec(end+1) = true;

% combine results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
