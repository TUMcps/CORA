function hs = groupHalfspaces(list,dom,varargin)
% groupHalfspaces - replace a list of similar halfspaces by a single
%                   halfspace that outer- or inner-approximates the
%                   intersection of the union of halfspaces with a domain
%
% Syntax:
%    hs = groupHalfspaces(list,dom)
%    hs = groupHalfspaces(list,dom,type)
%
% Inputs:
%    list - list of halfspaces (class halfspace) represented as cell-array
%    dom - domain of values (class contSet)
%    type - inner-approximation ("inner") or outer-approximation ("outer")
%
% Outputs:
%    hs - resulting halfspace
% 
% Example: 
%    dom = interval(-ones(2,1),ones(2,1));
%    hs1 = halfspace([-1 -1],-1);
%    hs2 = halfspace([-1.2 -0.9],-1);
%    hs3 = halfspace([-0.9 -1.1],-1);
%
%    hs_i = groupHalfspaces({hs1,hs2,hs3},dom,'inner');
%    hs_o = groupHalfspaces({hs1,hs2,hs3},dom,'outer');
%
%    figure; hold on; box on;
%    plot(dom);
%    xlim([-2,2]); ylim([-2,2]);
%    plot(hs1,[1,2],'r','FaceAlpha',0.5);
%    plot(hs2,[1,2],'r','FaceAlpha',0.5);
%    plot(hs3,[1,2],'r','FaceAlpha',0.5);
%    plot(conHyperplane(hs_i),[1,2],'g');
%    plot(conHyperplane(hs_o),[1,2],'c');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Authors:       Niklas Kochdumper
% Written:       09-November-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    % parse input arguments
    type = 'outer';
    
    if nargin > 2 && ~isempty(varargin{1})
        type = varargin{1};
    end

    % compute intervals for normal vector and offset
    n = norm(list{1}.c);
    c = interval(list{1}.c / n); d = interval(list{1}.d / n);

    for i = 2:length(list)
        n = norm(list{1}.c);
        c = c | (list{i}.c / n); d = d | (list{i}.d / n);
    end

    % construct new halfspace containing all halfspaces in the list
    if ~isa(dom,'interval')
        dom = interval(dom);
    end

    c_ = center(c); d_ = d - (c - c_)'*dom;

    if strcmp(type,'outer')
        hs = halfspace(c_,supremum(d_));
    else
        hs = halfspace(c_,infimum(d_));
    end
end

% ------------------------------ END OF CODE ------------------------------
