classdef testClass
    %TEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        a
    end
    
    methods
        function obj = testClass(a)
            obj.a = a;
        end
        
        function y = eval(obj,x)
           y = obj.a*x'*x; 
        end
        
        function g = diff(obj,x)
           g = 2*obj.a*x'; 
        end
    end
    
end

