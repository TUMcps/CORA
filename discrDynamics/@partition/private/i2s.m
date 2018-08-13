function subscriptMatrix=i2s(Obj,indexVector)
% Purpose:  convert segment number into segment subscripts
% Pre:      1st Parameter - partition object
%           2nd Parameter - row vector of cell indices; 
% Post:     Return segment subscripts
% Tested:   14.09.06,MA
% Checked:  1.8.17 AP

    %obtain vector of number of segments for each dimension as a row vector
    siz=Obj.nrOfSegments';

    %subscript variable string
    string=[];
    for iChar=1:length(siz)
        string=[string,'s',num2str(iChar),','];
    end
    string(end)=[];
    %Generate command string
    command=['[',string,']=ind2sub([',num2str(siz),'],[',num2str(indexVector),']);'];
    eval(command);

    %arrange variables in a vector
    string=strrep(string,',',';');
    command=['subscriptMatrix=[',string,']'';'];
    eval(command);
    end
%end