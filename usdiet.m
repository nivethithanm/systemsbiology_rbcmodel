%custom function for incorporating diets
function v = usdiet(model,diet)   
exc = 0;
for i=1:length(model.rxns)
    if isequal(model.rxns{i}(1:2),'EX')
        exc = exc + 1;
    end
end
for i=1:length(diet)
    name = diet{i,1};
    for k=1:length(name)
        if isequal(name(k),'[')
            name(k) = '(';
        end
        if isequal(name(k),']')
            name(k) = ')';
        end
    end
    for j=1:exc
        if isequal(name,model.rxns{j,1})
            model.lb(j) = -diet{i,2};
        end
    end
end
v = model;
end
