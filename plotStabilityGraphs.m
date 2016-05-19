addpath('continuity_netsci')
addpath('packages/mcode')
addpath('causaldisc2016/')

nodelabels = readtable('pcnets/roinames.csv','ReadVariableNames',1, 'Delimiter',',');
nodelabels = nodelabels.Properties.VariableNames;
dataPath = 'causaldisc2016/Data/figures/'


% Method:PC -------------------------------------------------------------------
method = 'PC'
basefilename = 'Resting_Stability_PC'
inputFile = 'causaldisc2016/Data/Resting_Stability_PC_2016_05_18';
testMatrix = readtable(inputFile, 'ReadRowNames',1,'ReadVariableNames',1);
% node of interest
figobj = figure(1);
set(figobj,'Position',[440 122 803 676]);
NOI = 'LAMFG'
try
	dir([dataPath  method]);
catch
	mkdir([dataPath  method]);
end
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);

figure(1)
% node of interest
NOI = 'LPMFG'
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);

figure(1)
% node of interest
NOI = 'RAMFG'
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);

figure(1)
% node of interest
NOI = 'RPMFG'
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);


% Method: GIES
method = 'GIES'
basefilename = 'RestingRight_Stability_GIES'
inputFile = ['causaldisc2016/Data/GIES/RightStim_Stability_GIES'];
dataPath = 'causaldisc2016/Data/figures/'
load(inputFile,'stabMat'); p = size(stabMat,2);stabMat = reshape(stabMat,[p p size(stabMat,1)/p]);
testMatrix = array2table(mean(stabMat,3));

figno = 2;
figobj = figure(figno);
set(figobj,'Position',[440 122 803 676]);
NOI = 'LAMFG'
try
	dir([dataPath  method]);
catch
	mkdir([dataPath  method]);
end
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);

figure(figno)
% node of interest
NOI = 'LPMFG'
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);

figure(figno)
% node of interest
NOI = 'RAMFG'
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);

figure(figno)
% node of interest
NOI = 'RPMFG'
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);



basefilename = 'RestingLeft_Stability_GIES'
inputFile = ['causaldisc2016/Data/GIES/LeftStim_Stability_GIES'];
dataPath = 'causaldisc2016/Data/figures/'
load(inputFile,'stabMat'); p = size(stabMat,2);stabMat = reshape(stabMat,[p p size(stabMat,1)/p]);
testMatrix = array2table(mean(stabMat,3));

figno = 2;
figobj = figure(figno);
set(figobj,'Position',[440 122 803 676]);
NOI = 'LAMFG'
try
	dir([dataPath  method]);
catch
	mkdir([dataPath  method]);
end
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);

figure(figno)
% node of interest
NOI = 'LPMFG'
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);

figure(figno)
% node of interest
NOI = 'RAMFG'
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);

figure(figno)
% node of interest
NOI = 'RPMFG'
[h] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,[dataPath method '/']);
