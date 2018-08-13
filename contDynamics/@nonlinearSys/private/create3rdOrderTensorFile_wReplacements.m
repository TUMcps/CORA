function create3rdOrderTensorFile_wReplacements(J3dyn,path,rep)
% create3rdOrderTensorFile - generates an mFile that allows to compute the
% 3rd order terms 
%
% Syntax:  
%    createHessianTensorFile(obj,path)
%
% Inputs:
%    obj - nonlinear system object
%    path - path for saving the file
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
% Written:      22-August-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


for k=1:length(J3dyn(:,1,1,1))
    for l=1:length(J3dyn(1,:,1,1))
        Tdyn{k,l} = squeeze(J3dyn(k,l,:,:));
    end
end

fid = fopen([path '/thirdOrderTensor.m'],'w');
fprintf(fid, '%s\n\n', 'function [T]=thirdOrderTensor(x,u)');


%write replacements
for iRep = 1:length(rep)
    strVar = char(rep{iRep}.var);
    strVal = bracketSubs(rep{iRep}.val);
    str = [strVar,' = ',strVal,';'];
    %write in file
    fprintf(fid, '%s\n', str);
end


%create list to compare equals
counter = 1;
for k=1:length(Tdyn(:,1))
    for l=1:length(Tdyn(1,:))
        M = Tdyn{k,l};
        for iRow = 1:(length(M(:,1)))
            for iCol = 1:(length(M(1,:)))
                list(counter).val = M(iRow,iCol);
                list(counter).ind = [k l iRow iCol];
                counter = counter + 1;
            end
        end
    end
end

%write part for parallel execution
dim = length(Tdyn(:,1));
str = ['C = cell(',num2str(dim),',1);'];
%write in file
fprintf(fid, '\n\n %s\n\n', str);
str = ['parfor i = 1:',num2str(dim)];
%write in file
fprintf(fid, '\n\n %s\n\n', str);
str = 'switch i';  
%write in file
fprintf(fid, '\n\n %s\n\n', str);



%dynamic part
notReplacedList = [];
for k=1:length(Tdyn(:,1))
    caseStr = ['case ',num2str(k)];
    %write in file
    fprintf(fid, '\n\n %s\n\n', caseStr);
    for l=1:length(Tdyn(1,:))
        %get matrix size
        [rows,cols] = size(Tdyn{k,l});
        sparseStr = ['sparse(',num2str(rows),',',num2str(cols),')'];
        %str=['T{',num2str(k),',',num2str(l),'} = infsup(',sparseStr,',',sparseStr,');'];
        str=['C{i}{',num2str(l),'} = interval(',sparseStr,',',sparseStr,');'];
        %write in file
        fprintf(fid, '\n\n %s\n\n', str);
        % write rest of matrix
        %notReplacedList = writeSparseMatrix(Tdyn{k,l},['T{',num2str(k),',',num2str(l),'}'],fid,rep,list,k,l,notReplacedList);
        notReplacedList = writeSparseMatrix(Tdyn{k,l},['C{',num2str(k),'}{',num2str(l),'}'],fid,rep,list,k,l,notReplacedList);

        disp(['dynamic index ',num2str(k),',',num2str(l)]);
    end
end


%write end of parallel execution
fprintf(fid, '\n\n %s\n', 'end');
fprintf(fid, '%s\n\n', 'end');
str = ['for i=1:',num2str(dim)];
fprintf(fid, '%s\n', str);
str = ['for j=1:',num2str(length(Tdyn(1,:)))];
fprintf(fid, '%s\n', str);
str = ['T{i,j} = C{i}{j};'];
fprintf(fid, '%s\n', str);
fprintf(fid, '%s\n', 'end');
fprintf(fid, '%s\n', 'end');



%close file
fclose(fid);


function [notReplacedList] = writeSparseMatrix(M,var,fid,rep,list,kInd,lInd,notReplacedList)


%write each row
for iRow = 1:(length(M(:,1)))
    for iCol = 1:(length(M(1,:)))
        if M(iRow,iCol) ~= 0
%             %replacements
%             for iRep = 1:length(rep)
%                 M(iRow,iCol) = subs(M(iRow,iCol), rep{iRep}.val, rep{iRep}.var);
%             end

            %initialize
            try
                str=char(M(iRow,iCol));
            catch
                str=sprintf('%f',M(iRow,iCol));
            end
            
            %check for equal values
            keepGoing = 1;
            iInd = 1;
            counter = 1;
            replaced = 0;
            while keepGoing
                if isequal(M(iRow,iCol), list(iInd).val)
                %if strcmp(char(M(iRow,iCol)), char(list(iInd).val))
                    keepGoing = 0;
                    if ~all(list(iInd).ind == [kInd, lInd, iRow, iCol])
                        %str=['T{',num2str(list(iInd).ind(1)),',',num2str(list(iInd).ind(2)),'}(',num2str(list(iInd).ind(3)),',',num2str(list(iInd).ind(4)),')'];
                        str=['C{i}{',num2str(list(iInd).ind(2)),'}(',num2str(list(iInd).ind(3)),',',num2str(list(iInd).ind(4)),')'];
                        replaced = 1;
                    else
                        notReplacedList(end+1) = iInd;
                    end
                end
                counter = counter + 1;
                try
                    iInd = notReplacedList(counter);
                catch
                    iInd = iInd + 1;
                end
            end

            %replacements
            if ~replaced
                %add spacing
                str = [str,' '];
                for iRep = 1:length(rep)
                    %remove with brackets
                    str = strrep(str, ['(',rep{iRep}.val,')'], char(rep{iRep}.var));
                    %remove stand alone
                    str = strrep(str, [' ',rep{iRep}.val,' '], char(rep{iRep}.var));
%                     %remove if only contains * and ^
%                     ind1 = find(rep{iRep}.val == '+');
%                     ind2 = find(rep{iRep}.val == '-');
%                     if isempty([ind1 ind2])
%                         %remove attached to a multiplication
%                         str = strrep(str, ['*',rep{iRep}.val], ['*',char(rep{iRep}.var)]);
%                     end
                end

                %remove spacing at end
                if str(end) == ' '
                    str(end)=[];
                end
            end
            
            %remove brackets
            str=bracketSubs(str);
            
            str=[var,'(',num2str(iRow),',',num2str(iCol),') = ',str,';'];
            %write in file
            fprintf(fid, '%s\n', str);
        end
    end
end



function [str]=bracketSubs(str)

%generate left and right brackets
str=strrep(str,'L','(');
str=strrep(str,'R',')');

%------------- END OF CODE --------------