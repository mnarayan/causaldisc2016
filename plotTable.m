nodelabels = readtable('pcnets/roinames.csv','ReadVariableNames',1, 'Delimiter',',');
nodelabels = nodelabels.Properties.VariableNames;
dataPath = 'causaldisc2016/Data/'

% Method: GIES
method = 'GIES';
method_num = 1;
basefilename = 'RestingRight_Stability_GIES'
inputFile = ['causaldisc2016/Data/GIES/RightStim_Stability_GIES'];
dataPath = 'causaldisc2016/Data/figures/'
load(inputFile,'stabMat'); p = size(stabMat,2);stabMat = reshape(stabMat,[p p size(stabMat,1)/p]);
testMatrix = array2table(mean(stabMat,3));
testMatrix.Properties.RowNames = nodelabels; 
testMatrix.Properties.VariableNames = nodelabels;

selected_NOIS = {'LAMFG','LPMFG','RAMFG','RPMFG'};
tblresults1 = table(); tblresults2 = table();
tbl_rownames = {};
counter = 0;
for ii=1:length(selected_NOIS)
	for jj=ii:length(selected_NOIS)
		if(ii~=jj)
			counter=counter+1;
			tblresults1(counter,method_num) = testMatrix(selected_NOIS(ii),selected_NOIS(jj));
			tbl_rownames = cat(2, tbl_rownames, {strcat(selected_NOIS(ii), '_2_' ,selected_NOIS(jj))});
		end
	end
end
counter = 0;
for ii=1:length(selected_NOIS)
	for jj=ii:length(selected_NOIS)
		if(ii~=jj)
			counter = counter+1;
			tblresults2(counter,method_num) = testMatrix(selected_NOIS(jj),selected_NOIS(ii));
			tbl_rownames = cat(2, tbl_rownames, {strcat(selected_NOIS(jj), '_2_' ,selected_NOIS(ii))});
		end
	end
end
tblresults = vertcat(tblresults1,tblresults2);
tblresults.Properties.VariableNames = {basefilename};
tblresults.Properties.RowNames = [tbl_rownames{:}];


%%%%%%%%%%%%%%
method_num = 2;
basefilename = 'RestingLeft_Stability_GIES'
inputFile = ['causaldisc2016/Data/GIES/LeftStim_Stability_GIES'];
dataPath = 'causaldisc2016/Data/figures/'
load(inputFile,'stabMat'); p = size(stabMat,2);stabMat = reshape(stabMat,[p p size(stabMat,1)/p]);
testMatrix = array2table(mean(stabMat,3));
testMatrix.Properties.RowNames = nodelabels; 
testMatrix.Properties.VariableNames = nodelabels;

selected_NOIS = {'LAMFG','LPMFG','RAMFG','RPMFG'};
tblresults1 = table(); tblresults2 = table();
tbl_rownames = {};
counter = 0;
for ii=1:length(selected_NOIS)
	for jj=ii:length(selected_NOIS)
		if(ii~=jj)
			counter=counter+1;
			tblresults1(counter,1) = testMatrix(selected_NOIS(ii),selected_NOIS(jj));
			tbl_rownames = cat(2, tbl_rownames, {strcat(selected_NOIS(ii), '_2_' ,selected_NOIS(jj))});
		end
	end
end
counter = 0;
for ii=1:length(selected_NOIS)
	for jj=ii:length(selected_NOIS)
		if(ii~=jj)
			counter = counter+1;
			tblresults2(counter,1) = testMatrix(selected_NOIS(jj),selected_NOIS(ii));
			tbl_rownames = cat(2, tbl_rownames, {strcat(selected_NOIS(jj), '_2_' ,selected_NOIS(ii))});
		end
	end
end
tblresults = horzcat(tblresults,vertcat(tblresults1,tblresults2));
tblresults.Properties.VariableNames(method_num) = {basefilename};
[tbl_rownames{:}];
writetable(tblresults,[dataPath  method '_table' datestr(now,'dd.mm.yyyy') '.csv'],'Delimiter',',');

% Now use this table as input in our input struct:
input.data = tblresults;
% Set the row format of the data values (in this example we want to use
% integers only):
input.dataFormat = {'%.3f'};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.tableBorders = 1;
% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;
% Now call the function to generate LaTex code:
latex = latexTable(input);
% write to file
fileID = fopen([dataPath  method '_table' datestr(now,'dd.mm.yyyy') '_table.tex'],'w');
formatSpec = '%s\n';
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fileID,formatSpec,latex{row,:});
end
fclose(fileID);
type([dataPath  method '_table' datestr(now,'dd.mm.yyyy') '_table.tex'])
