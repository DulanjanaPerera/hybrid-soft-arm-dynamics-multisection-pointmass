function state=residualPressureInitialState()
%#codegen
state=pacPressureInitialState();state.residualAnchor=zeros(2,1);
state.linearRequest=zeros(2,1);state.lastEffective=zeros(2,1);
end
