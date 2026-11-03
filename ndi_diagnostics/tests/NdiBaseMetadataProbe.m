classdef NdiBaseMetadataProbe < AuroraReader_s2_BaseRotation
    methods
        function [sizes,types,fixed,complex]=metadata(obj)
            sizes=cell(1,23);[sizes{:}]=getOutputSizeImpl(obj);
            types=cell(1,23);[types{:}]=getOutputDataTypeImpl(obj);
            fixed=cell(1,23);[fixed{:}]=isOutputFixedSizeImpl(obj);
            complex=cell(1,23);[complex{:}]=isOutputComplexImpl(obj);
        end
    end
end
