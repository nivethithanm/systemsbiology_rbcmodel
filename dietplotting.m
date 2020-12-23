rbcmodel1 = readCbModel('RBC.xml');

testfiledir = 'Diets\';
matfiles = dir(fullfile(testfiledir, '*.xlsx'));
nfiles = length(matfiles);
afcrxns = {};
Jc = [];
for i = 1 : nfiles
dietConstraints = readtable(['Diets/' matfiles(i).name]);
dietConstraints=table2cell(dietConstraints);
rbcmodel1 = usdiet(rbcmodel1,dietConstraints);
[min0,max0] = fluxVariability(rbcmodel1,0);

% genes list = {Gpi.1, {{'Ldha.1';'Ldhb.1'}, 'Taldo1.1'}}
gene = 'No Gene Deletion';

rbcdel1 = deleteModelGenes(rbcmodel1, 'Taldo1.1');  

[min1,max1] = fluxVariability(rbcdel1,0);

J = fvaJaccardIndex([min0, min1], [max0, max1]);
Jc = [Jc J];

E = [(max0 - min0)/2 (max1 - min1)/2];
Y = [min0 min1] + E;
X = [(1:length(Y)) - 0.1; (1:length(Y)) + 0.1]';

%[~, xj] = sort(J);
xj = 1:length(J);
f1 = figure;
if strcmp(version('-release'), '2016b')
    errorbar(X, Y(xj, :), E(xj, :), 'linestyle', 'none', 'linewidth', 2, 'capsize', 0);
else
    %errorbar(X, Y(xj, :), E(xj, :), 'linestyle', 'none', 'linewidth', 2);
    hold on
    errorbar(X(:,1), Y(xj, 1), E(xj, 1), 'linestyle', 'none', 'linewidth', 2,'Color','b');
    errorbar(X(:,2), Y(xj, 2), E(xj, 2), 'linestyle', 'none', 'linewidth', 2,'Color','r');
end
set(gca, 'xlim', [0, length(Y) + 1])
xlabel('Reaction')
ylabel('Flux range (mmol/gDW/h)')
ylim([-50,50])
yyaxis right
%plot(J(xj),'linewidth', 2)
legend('Normal', 'Deleted', 'Jaccard','location', 'northoutside', ...
       'orientation', 'horizontal')
title({['diet : ', matfiles(i).name(1:(length(matfiles(i).name)-5))],['gene : ', gene]});
saveas(f1,['diet_', matfiles(i).name(1:(length(matfiles(i).name)-5)), gene,'.jpg']);
%}
afcrxns{i,1} = threshold(rbcmodel1,rbcdel1,0.7);   %threshold similarity = 70%
afcrxns{i,2} = matfiles(i).name(1:(length(matfiles(i).name)-5));
fprintf('progres %f',i);
end
surfJc);  %3D plot