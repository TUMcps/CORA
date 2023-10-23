function simRes = add(simRes1,simRes2)
% add - joins two simResult objects
%
% Syntax:
%    simRes = add(simRes1,simRes2)
%
% Inputs:
%    simRes1 - simResult object
%    simRes2 - simResult object
%
% Outputs:
%    simRes - resulting simResult object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Niklas Kochdumper
% Written:       29-May-2020             
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% special case
if isempty(simRes1)
     simRes = simRes2;
     return;
elseif isempty(simRes2)
     simRes = simRes1;
     return;
end

% check input arguments
inputArgsCheck({{simRes1,'att',{'simResult'}};
                {simRes2,'att',{'simResult'}}});

% general case
if isempty(simRes1.loc) ~= isempty(simRes2.loc)
    throw(CORAerror('CORA:specialError','Objects are not compatible!')); 
end

simRes = simRes1;
simRes.x = [simRes1.x; simRes2.x];
simRes.t = [simRes1.t; simRes2.t];
simRes.loc = [simRes1.loc; simRes2.loc];

% ------------------------------ END OF CODE ------------------------------
