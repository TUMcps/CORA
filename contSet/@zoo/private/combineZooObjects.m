function [obj1,obj2] = combineZooObjects(obj1,obj2)


    % find all methods that appear in both zoo-objects
    objects1 = cell(length(obj1.objects),1);
    objects2 = cell(length(obj2.objects),1);
    methods = cell(length(obj1.objects),1);
    
    counter = 1;
    
    for i = 1:length(obj1.method)
        [~,ind] = ismember(obj1.method{i},obj2.method);
        if ind > 0
            objects1{counter} = obj1.objects{i};
            objects2{counter} = obj2.objects{ind};
            methods{counter} = obj1.method{i};
            counter = counter + 1;
        end
    end

    % construct the adapted zoo-objects
    obj1.method = methods(1:counter-1);
    obj1.objects = objects1(1:counter-1);
    obj2.objects = objects2(1:counter-1);