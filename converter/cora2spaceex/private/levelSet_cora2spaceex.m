function eqs = levelSet_cora2spaceex(obj)
% levelSet_cora2spaceex - generate a string that decribes the object
%                         in SpaceEx format
%
% Syntax:
%    eqs = levelSet_cora2spaceex(obj)
%
% Inputs:
%    obj - levelSet object
%
% Outputs:
%    eqs - string describing the set in SpaceEx format
%
% Example:
%    syms x y
%    eq = x^2 + y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%
%    levelSet_cora2spaceex(ls)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: cora2spaceex

% Authors:       Niklas Kochdumper
% Written:       18-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    eqs = [];

    % get object properties
    compOp = obj.compOp;
    if ~iscell(compOp)
       compOp = {compOp}; 
    end
    
    x = sym('x',[length(obj.vars),1]);
    eq = obj.funHan(x);
    
    % loop over all equations
    for i = 1:length(compOp)
        
        % different comparison operators
        if strcmp(compOp{i},'<=')
            temp = [char(eq(i)),' <= 0'];
        elseif strcmp(compOp{i},'<')
            temp = [char(eq(i)),' < 0'];
        else
            temp = [char(eq(i)),' == 0'];
        end
        
        % add current equation to overall string
        if isempty(eqs)
            eqs = temp; 
        else
            eqs = [eqs,' & ',newline,temp];
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
