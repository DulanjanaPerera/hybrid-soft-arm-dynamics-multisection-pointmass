function p=sgPlain(p)
for f=fieldnames(p)'
    if isa(p.(f{1}),'dlarray'),p.(f{1})=extractdata(p.(f{1}));end
end
end
