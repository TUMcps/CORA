function nameMat = addMatIndizesToName(name,mat)
    
   x = 1:size(mat,1);
   y = 1:size(mat,2);
   [X,Y] = meshgrid(x,y);
   Z = 10*X'+Y';

   temp1 = cellfun(@(x) num2str(x),num2cell(Z),'UniformOutput',false);
   temp2 = repmat({name},size(mat));
   nameMat = cellfun(@(x,y) strcat(x,y),temp2,temp1,'UniformOutput',false);



end