rbcmodel1 = readCbModel('RBC.xml');

[min0,max0] = fluxVariability(rbcmodel1,0);

% genes list = {Gpi.1, {{'Ldha.1';'Ldhb.1'}, 'Taldo1.1'}}
gene = 'Taldo1.1';

rbcdel1 = deleteModelGenes(rbcmodel1,'Taldo1.1');  

[min1,max1] = fluxVariability(rbcdel1,0);

J = fvaJaccardIndex([min0, min1], [max0, max1]);
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
%ylabel('Jaccard index')
title({['diet : ','No Diet'],['gene : ', gene]});
saveas(f1,['no_diet_', gene,'.jpg']);