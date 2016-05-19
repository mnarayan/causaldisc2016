function [h2] = plot_schemadigraph(testMatrix,NOI,nodelabels,basefilename,dataPath)


	% Extract local neighborhood
	NOI_idx = ~cellfun(@isempty, strfind(nodelabels,NOI));
	node_idx = find(NOI_idx); neighbor_idx = find(~NOI_idx);
	nrois = size(testMatrix,1);
	subGraph = zeros(size(table2array(testMatrix)));
	subGraph(node_idx,neighbor_idx) = table2array(testMatrix(node_idx,neighbor_idx));
	subGraph(neighbor_idx,node_idx) = table2array(testMatrix(neighbor_idx,node_idx));
	
	left_regions = [1:12]; 
	reverse_idx2 = fliplr([left_regions]); % left/right bilateral symmetry;
	tmp_subGraph = subGraph;
	tmp_subGraph(left_regions,left_regions) = subGraph(reverse_idx2,reverse_idx2); 
	tmp_subGraph(left_regions,setdiff(1:30,left_regions)) = subGraph(reverse_idx2,setdiff(1:30,left_regions));
	tmp_subGraph(setdiff(1:30,left_regions),left_regions) = subGraph(setdiff(1:30,left_regions),reverse_idx2);
	subGraph = tmp_subGraph; 
	nodelabels = [nodelabels(reverse_idx2) nodelabels(setdiff(1:30,left_regions))];
	clear tmp_subGraph;
	
	tmp_undG = digraph2bipartite((subGraph));
	undG = tmp_undG; 
	reverse_idx = fliplr([31:60]);
	undG(1:30,31:60) = tmp_undG(1:30,reverse_idx);
	undG(31:60,1:30) = tmp_undG(reverse_idx,1:30);
	undG(31:60,31:60) = tmp_undG(reverse_idx,reverse_idx);
	clear tmp_undG;

	% Orange vs Blue
	% NodeColor =  cat(1,repmat([0 60 50]/100, [30 1]), repmat([80 50 10]/100, [30 1]));
	% Light vs Dark Orange
	%NodeColor =  cat(1,repmat([80 50 10]/100, [30 1]), repmat([40 25 10]/100, [30 1]));
	% Light vs Dark Blue
	NodeColor =  cat(1,repmat([0 60 50]/100, [30 1]), repmat([0 40 30]/100, [30 1]));
	maincolor = [.6 .2 .8]; 
	fontsize = 16;
	
%	useSchemanet = 1;
%	if(useSchemanet)
	h2 = schemanet(undG,[nodelabels fliplr(nodelabels)],[0 1],'copper',2,{[1:30]'; [31:60]';}); set(gcf,'Color',[.9 .9 .9]);
	set(h2.cb,'Color',[.1 .1 .1]);
	set(h2.o(1),'Color',NodeColor(1,:)); set(h2.o(2),'Color',NodeColor(31,:));
	if(strfind(version,'2015')|strfind(version,'2016'))
		
		set(h2.s,'visible','off');
		node_x = get(h2.s,'XData'); node_y = get(h2.s,'YData'); 
		h2.s(1) = scatter(node_x(1:30),node_y(1:30),'o','SizeData',75,'MarkerEdgeColor',[.1 .1 .1],'LineWidth',2);
		h2.s(2) = scatter(node_x(31:60),node_y(31:60),'x','SizeData',100,'MarkerEdgeColor',[.1 .1 .1],'LineWidth',2);
		set(h2.s(1), 'CData',repmat([.1 .1 .1],[60 1]));
		set(h2.s(1), 'CData', NodeColor(1:30,:));
		set(h2.s(2), 'CData',repmat([.1 .1 .1],[60 1]));
		set(h2.s(2), 'CData', NodeColor(31:60,:));
		
		%set(h2.s,'Marker',{cat(1,repmat('o',[30 1]),repmat('x',[30 1]))});
        %set(h2.s,'MarkerEdgeColor',[.3 .3. .3],'MarkerSize',15,'Marker',cat(1,repmat('o',[30 1]),repmat('x',[30 1])));		
		t_idx = find(~isnan(h2.t));		
		for ii=1:length(t_idx)
			ii;			
			set(h2.t(t_idx(ii)),'Color',NodeColor(ii,:),'FontSize',fontsize);
			% if(ii<=30)
			%   set(h2.s, 'CData',NodeColor','Marker','x', 'MarkerSize',15,'MarkerEdgeColor',[.3 .3 .3],'MarkerFaceColor',[.3 .3 .3]);
			% else
			%   set(h2.s, 'CData',NodeColor','Marker', 'o', 'MarkerSize',12,'MarkerEdgeColor',[.3 .3 .3],'MarkerFaceColor',[.3 .3 .3]);
			% end
		end
	else
		h2.n = get(h2.s,'Children');
		t_idx = find(~isnan(h2.t));		
		for ii=1:length(t_idx)
			ii;
			set(h2.t(t_idx(ii)),'Color',NodeColor(ii,:));
			if(ii<=30)
			  set(h2.n(t_idx(ii)), 'Marker','x', 'MarkerSize',15,'MarkerEdgeColor',[.3 .3 .3],'MarkerFaceColor',[.3 .3 .3]);
			else
			  set(h2.n(t_idx(ii)), 'Marker', 'o', 'MarkerSize',12,'MarkerEdgeColor',[.3 .3 .3],'MarkerFaceColor',[.3 .3 .3]);
			end
		end
	end
	l_idx = find(~isnan(h2.l));
	max_lw = 7;
	for jj =1:length(l_idx)
		tmp_col = get(h2.l(l_idx(jj)),'Color');
		set(h2.l(l_idx(jj)),'LineWidth', max(1,max_lw*tmp_col(1)));
	end
	addpath(genpath('packages/matlab2tikz'))
	set(h2.o,'Visible','off');
	set(gcf,'PaperPosition',[0.2500 2.5000 8 8],'PaperOrientation','Landscape');
	savefig([dataPath 'Digraph_' NOI '_' basefilename datestr(now,'dd.mm.yyyy') '.fig']);
	print('-dpng','-r600', [dataPath 'Digraph_' NOI '_' basefilename datestr(now,'dd.mm.yyyy') '.png'])
	matlab2tikz('filename',[dataPath 'Digraph_' NOI '_' basefilename '.tikz'],'floatFormat', '%.3f','externalData', false, ...
	                    'height', '.3\textwidth','width', '.3\textwidth');
% 	  else
% 	  h = schemaball(undG, ...
% 	          [nodelabels fliplr(nodelabels)], ...
% 	  			  maincolor, [1 1 1]);
% 	          %strcat( vertcat(repmat({'L.','R.'}',[28 1]),repmat({'L.'},[7 1]),repmat({'R.'},[7 1]),repmat({'L.','R.'}',[20 1])),SortedROIs), ...
% 	  set(h.s, 'MarkerEdgeColor','black','LineWidth',3,'SizeData',125,'CData',ones(length(undG),1)*[.3 .3 .3])
% 	  h.n = get(h.s,'Children')
% 	  %set(h.l(~isnan(h.l)), 'LineWidth',3);
% 	  set(h.l(setdiff(find(~isnan(h.l)), [1])), 'LineWidth',3, 'Color', maincolor);
% 	  set(gca,'Color', [1 1 1]); axis off;
% 	  set(h.l(1),'Color',[.7 .7 .7])
% 	  set(h.t,'Color',[.2 .2 .2],'fontsize',20,'FontName','Helvetica','FontWeight','bold');
% 	  t_idx = find(~isnan(h.t));
% 	  for ii=1:length(t_idx)
% 	  	ii
% 	    set(h.t(t_idx(ii)),'Color',NodeColor(ii,:));
% 	    if(ii>30)
% 	      set(h.n(ii), 'Marker','o');
% 	    else
% 	      set(h.n(ii), 'Marker', 'x');
% 	    end
% 	  end
% 	end
end