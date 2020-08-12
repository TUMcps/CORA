function options = set_inputSet(obj, options, checkName)
% set_inputSet - sets the interval options "uTrans", "uTransVec" and
%                "originContained"
%
% Syntax:
%    options = set_inputSet(obj, options, checkName)
%
% Inputs:
%    obj       - contDyn object
%    options   - options for object
%    checkName - checkOptions* function (for error tracing)
%
% Outputs:
%    options   - options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: check_inputSet

% Author:       Niklas Kochdumper, Matthias Althoff
% Written:      02-January-2020
% Last update:  07-January-2020 (including constant input)
%               17-July-2020 (linProbSys added, MA)
% Last revision:---

%------------- BEGIN CODE --------------

    % set options "originContained"
    if ~strcmp(checkName,'checkOptionsSimulate') && ...
       (isa(obj,'linearSys') || isa(obj,'linParamSys') || isa(obj,'linProbSys'))

        if isfield(options,'u')
           options.originContained = 0;
        elseif isa(obj,'linearSys') && ~isempty(obj.c)
            % check if any of equality constraints B(i)*U = -c(i) = 0
            % is unsatisfiable
            found = 0;
            for i = 1:size(obj.B,1)
               hp = conHyperplane(obj.B(i,:),-obj.c(i));
               if ~isIntersecting(hp,options.U)
                  options.originContained = 0;
                  found = 1;
               end
            end
            if ~found
                % check if B*U + c (=vTrans) contains 0         
                vTrans = obj.B*options.U + obj.c;
                % faster computation if vTrans is an interval
                if isInterval(vTrans)
                    vTransInt = interval(vTrans);
                    options.originContained = in(vTransInt,zeros(dim(vTrans),1));
                else
                    options.originContained = in(vTrans,zeros(dim(vTrans),1));
                end
            end
        elseif isInterval(options.U)
            % no constant input, faster computation if U is an interval
            int = interval(options.U);
            options.originContained = in(int,zeros(dim(options.U),1));
        else
            % no constant input c, U not an interval
            options.originContained = in(options.U,zeros(dim(options.U),1));
        end
    end

    % set options "uTrans" and "uTransVec"
    uTrans = center(options.U);
    options.U = options.U - uTrans;

    if isfield(options,'u')
        options.uTransVec = options.u + uTrans;
        options = rmfield(options,'u');
    else
        options.uTrans = uTrans;
    end
    
    
end

%------------- END OF CODE --------------