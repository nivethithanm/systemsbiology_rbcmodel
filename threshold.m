%finding affected reactions
function v = threshold(model0,model1,thsld)
[min0, max0] = fluxVariability(model0,0);
[min1, max1] = fluxVariability(model1,0);
J = fvaJaccardIndex([min0, min1], [max0, max1]);
%{
for i=0:0.001:1
    afcrxns = [];
    for j=1:76
        if J(j) < i
            afcrxns = [afcrxns j];
        end
    end
    if length(afcrxns) == 53
        fprintf('threshold= %f',i);
        break
    end
end
%}
afcrxns = [];
n = length(model0.rxns);
for j=1:n
    if J(j) < thsld
        afcrxns = [afcrxns model0.rxns(j)];
    end
end
v = afcrxns;
end