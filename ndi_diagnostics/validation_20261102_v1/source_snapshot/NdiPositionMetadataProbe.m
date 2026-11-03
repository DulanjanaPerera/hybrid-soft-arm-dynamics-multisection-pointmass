classdef NdiPositionMetadataProbe < AuroraReader_s2_Position_only
    methods
        function [sizes,types,fixed,complex]=metadata(obj)
            sizes=cell(1,20);[sizes{:}]=getOutputSizeImpl(obj);
            types=cell(1,20);[types{:}]=getOutputDataTypeImpl(obj);
            fixed=cell(1,20);[fixed{:}]=isOutputFixedSizeImpl(obj);
            complex=cell(1,20);[complex{:}]=isOutputComplexImpl(obj);
        end
    end
end
