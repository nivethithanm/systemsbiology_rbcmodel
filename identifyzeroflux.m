rbcmodel1 = readCbModel('RBC.xml');

testfiledir = 'Diets\';
matfiles = dir(fullfile(testfiledir, '*.xlsx'));
nfiles = length(matfiles);
afcrxns = {};
uzF = {};        % finding reactions with zero flux span
increaseJ = {};  % finding affected reactions in ascending J for finding biomarkers
for i = 1 : nfiles
dietConstraints = readtable(['Diets/' matfiles(i).name]);
%dietConstraints = readtable('Diets/DACH.xlsx');
dietConstraints=table2cell(dietConstraints);
rbcmodel1 = usdiet(rbcmodel1,dietConstraints);
[min0,max0] = fluxVariability(rbcmodel1,0);

% genes list = {'Gpi.1', {'Ldha.1';'Ldhb.1'}, 'Taldo1.1'}}
gene = 'No Gene Deletion';

rbcdel1 = deleteModelGenes(rbcmodel1,'Taldo1.1');  

%rbcdel1=usdiet(rbcmodel1,dietConstraints);

[min1,max1] = fluxVariability(rbcdel1,0);
%gene names = 'AmpD3.1' Bpgm.1 Bpgm.2 Rhag.1-e
uzf = [];             %unique zero flux
for j=1:length(min1)
    if min1(j) == 0
        if max1(j) == 0
            if ~((min0(j) == 0) & (max0(j) == 0))
                uzf = [uzf rbcmodel1.rxns(j)];
            end
        end
    end
end
uzF{i} = uzf;
J = fvaJaccardIndex([min0, min1], [max0, max1]);
[a,b] = sort(J);
c = [];
for j=1:99
c = [c; [rbcmodel1.rxns(j),J(j)]];
        
end
increaseJ{i} = c;
fprintf('progres %f',i);
end
%imagesc(Jc(1:(size(Jc,1))/2,:));