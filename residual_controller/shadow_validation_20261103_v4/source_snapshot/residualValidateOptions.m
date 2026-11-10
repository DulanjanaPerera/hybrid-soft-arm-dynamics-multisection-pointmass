function residualValidateOptions(options)
%#codegen
assert(all(isfinite(options))&&(options(1)==64||options(1)==128));
assert(options(2)>=1&&options(2)<=12&&options(2)==floor(options(2))&&options(3)>0);
assert(options(4)>=0&&options(5)>0&&all(options(6:8)>0));
end
