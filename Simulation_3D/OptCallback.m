classdef OptCallback<handle
    properties
        newiter = false
    end
    methods
        function stop = outfun(obj,x,optimValues,state)
            stop = false;
            if strcmp(state,'iter')
                obj.newiter = true;
            end
        end
    end
end