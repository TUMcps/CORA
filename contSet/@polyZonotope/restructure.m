function pZ = restructure(pZ,method,order,varargin)
% restructure - Calculates a new over-approxmiating representation of a
%    polynomial zonotope so that there remain no independent generators
%
% Syntax:
%    pZ = restructure(pZ, method, order)
%    pZ = restructure(pZ, method, order, genOrder)
%
% Inputs:
%    pZ - polyZonotope object
%    method - method used to calculate the new representation ('zonotope'
%             or 'reduce', 'reduceFull', or 'reducePart')
%    order - desired zonotope order of the dependent factors for the
%            resulting polynomial zonotope 
%    genOrder - desired zonotope order of the resulting polynomial zonotope
%               (only for method = 'reduce...')
%
% Outputs:
%    pZ - polyZonotope object over-approximating input polynomial zonotope
%
% Example:
%    pZ = polyZonotope([0;0],[1 0 1;1 2 -2],[-1 0.1 -0.5;1.2 0.3 0.2],[1 0 1;0 1 2]);
%    pZnew1 = restructure(pZ,'zonotopeGirard',2);
%    pZnew2 = restructure(pZ,'reduceGirard',2);
%
%    figure; hold on; ylim([-10,10]);
%    plot(pZnew2,[1,2],'g');
%    plot(pZ,[1,2],'r');
%
%    figure; hold on; ylim([-10,10]);
%    plot(pZnew1,[1,2],'b');
%    plot(pZ,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reduce

% Authors:       Niklas Kochdumper
% Written:       25-July-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
genOrder = setDefaultValues({Inf},varargin);

% check input arguments
inputArgsCheck({{pZ,'att','polyZonotope'};
                {order,'att','numeric','nonnan'};
                {method,'str',getMembers('restructureTechnique')};
                {genOrder,'att','numeric','nonempty'}});

% parse string for the method
if startsWith(method,'zonotope')
    spec = 1;
    redMeth = method(9:end);
elseif startsWith(method,'reduceFull')
    spec = 2;
    redMeth = method(11:end);
elseif startsWith(method,'reducePart')
    spec = 3;
    redMeth = method(11:end);
elseif startsWith(method,'reduceDI')
    spec = 4;
    redMeth = method(9:end);
elseif startsWith(method,'reduce')
    spec = 0;
    redMeth = method(7:end);
else
    throw(CORAerror('CORA:wrongValue','second',...
        "be 'zonotope', 'reduceFull', 'reducePart', 'reduceDI', or 'reduce'"));
end

redMeth(1) = lower(redMeth(1));

% restructure with the selected method
if spec == 1
    pZ = restructureZono(pZ,order,redMeth);
elseif spec == 2
    pZ = restructureReduceFull(pZ,order,redMeth);
elseif spec == 3
    pZ = restructureReducePart(pZ,order,redMeth);
elseif spec == 4
    pZ = restructureReduceDI(pZ,order,redMeth,genOrder);
else
    pZ = restructureReduce(pZ,order,redMeth,genOrder);
end

% ------------------------------ END OF CODE ------------------------------
