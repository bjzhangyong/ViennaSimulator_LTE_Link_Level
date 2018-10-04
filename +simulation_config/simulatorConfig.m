classdef simulatorConfig
    % This class defines a configuration class, which takes as input the
    % configuration struct of the simulator and adds certain parameters to
    % it.
    
    properties
    end
    
    methods (Static,Abstract)
        LTE_params = apply_parameters(LTE_params)
    end
end

