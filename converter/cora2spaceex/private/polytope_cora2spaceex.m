function eqs = polytope_cora2spaceex(obj)
% polytope_cora2spaceex - generate a string that decribes the object
%                            in SpaceEx format
%
% Syntax:
%    eqs = polytope_cora2spaceex(obj)
%
% Inputs:
%    obj - polytope object
%
% Outputs:
%    eqs - string describing the set in SpaceEx format
%
% Example:
%    A = [1 2;3 4;5 6];
%    b = [1;1;2];
%    poly = polytope(A,b);
%
%    polytope_cora2spaceex(poly)
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
        
    % inequality constraints A*x <= b
    A = obj.A;
    b = obj.b;

    x = sym('x',[size(A,2),1]);

    if ~isempty(A)

        Ax = A*x;

        for idx = 1: size(A,1)
            Astr = char(Ax(idx));
            bstr = num2str(b(idx));
            eq_c = [Astr, ' <= ' , bstr];
            if ~isempty(eqs)
                eqs = [eqs,' & ',newline,eq_c];
            else
                eqs = eq_c; 
            end
        end
    end

    % equality constraints Aeq*x = beq
    A = obj.Ae;
    b = obj.be;

    if ~isempty(A)

        Ax = A*x;

        for idx = 1: size(A,1)
            eq = [char(Ax(idx)), ' == ' , num2str(b(idx))];
            if ~isempty(eqs)
                eqs = [eqs,' & ',newline,eq];
            else
                eqs = eq; 
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
