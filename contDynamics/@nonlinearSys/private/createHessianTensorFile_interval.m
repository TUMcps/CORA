function createHessianTensorFile_interval(obj)
% createHessianTensorFile - generates an mFile that allows to compute the
% hessian tensor
%
% Syntax:  
%    createHessianTensorFile(obj)
%
% Inputs:
%    obj - nonlinear system object
%
% Outputs:
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      30-May-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    
%load second order jacobian
J2=obj.jacobians.secondOrder;

for k=1:length(J2(:,1,1))
    H{k} = squeeze(J2(k,:,:));
end


fid = fopen([coraroot '/contDynamics/stateSpaceModels/hessianTensor_interval.m'],'w');
fprintf(fid, '%s\n\n', 'function H=hessianTensor_interval(x,u)');
for k=1:length(H)
    str=['H{',num2str(k),'}=[...'];
    %write in file
    fprintf(fid, '%s\n', str);
    %save current matrix
    H_current = H{k};
    
    %loop through row vectors
    for iRow=1:length(H_current(1,:))
        str = char(H_current(iRow,1));
        %generate left and right brackets
        str = strrep(str,'L','(');
        str = strrep(str,'R',')');
        %replace zeros by interval-zeros
        if H_current(iRow,1)==0
            str=strrep(str,'0','interval(0,0)');
        end
        
        for iCol=2:length(H_current(:,1))
            %generate string
            strPartial = char(H_current(iRow,iCol));
            %generate left and right brackets
            strPartial = strrep(strPartial,'L','(');
            strPartial = strrep(strPartial,'R',')');
            %replace zeros by interval-zeros
            if H_current(iRow,iCol)==0
                strPartial=strrep(strPartial,'0','interval(0,0)');
            end
            %add to string
            str=[str,',',strPartial];
        end
        
        if iRow<length(H_current(1,:))
            %finalize string
            str=[str,'; ...'];
            %write in file
            fprintf(fid, '%s\n', str);
        else
            %finalize string
            str=[str,'];'];
            %write in file
            fprintf(fid, '%s\n', str);
        end
    end
end

%close file
fclose(fid);

%------------- END OF CODE --------------