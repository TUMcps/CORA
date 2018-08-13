function create3rdOrderTensorFileTaylorModel(J3dyn,path,name,x,u)
% create3rdOrderTensorFile - generates an mFile that allows to compute the
% 3rd order terms with taylor models
%
% Syntax:  
%    createHessianTensorFile(obj,path)
%
% Inputs:
%    J3dyn - symbolic third-order tensor
%    path - path for saving the file
%    name - name of the nonlinear function to which the 3rd order tensor belongs
%    x - symbolic variables for the states
%    u - symbolic variables for the inputs
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

% Author:       Niklas Kochdumper
% Written:      15-July-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


for k=1:length(J3dyn(:,1,1,1))
    for l=1:length(J3dyn(1,:,1,1))
        Tdyn{k,l} = squeeze(J3dyn(k,l,:,:));
    end
end

fid = fopen([path '/thirdOrderTensor_',name,'.m'],'w');
fprintf(fid, '%s\n\n', ['function [T]=thirdOrderTensor_',name,'(x,u,xInt,uInt)']);

%dynamic part
for k=1:length(Tdyn(:,1))
    for l=1:length(Tdyn(1,:))
        %get matrix size
        [rows,cols] = size(Tdyn{k,l});
        sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
        str=['T{',num2str(k),',',num2str(l),'} = interval(',sparseStr,',',sparseStr,');'];
        %write in file
        fprintf(fid, '\n\n %s\n\n', str);
        % write rest of matrix
        writeSparseMatrixTaylorModel(Tdyn{k,l},['T{',num2str(k),',',num2str(l),'}'],fid,x,u);

        disp(['dynamic index ',num2str(k),',',num2str(l)]);
    end
end

%close file
fclose(fid);

%------------- END OF CODE --------------