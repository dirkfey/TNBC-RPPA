%bmra_analysis_heatmaps ... Script to generate network heatmaps from
%   BMRAresults: 
%   1. NW heatmaps for all cell lines
%   2. Calculate consensus network for sensitive and resistant cells
%
%   Inputs: BMRAresults from bmra_run.m or saved in  mat_BMRAresults.mat
%   Outputs: 
%   1. '../figs/BMRI_NW_[CELLLINE]_3.png'
%   2. '../figs/BMRI_consNW_Sensitive.png', '../../figures/BMRI_consNW_Resistant.png'
%   Dependencies: df_pcolor, bmra_prepare_data

%% Load data (or alternative run the pipeline: bmra_run.m)
bmra_prepare_data           %get the network nodes, labels etc
datafile = 'mat_BMRAresults_19-09-09.mat';
load(datafile)    %load the results from bmra_run.m


%% Network heatmaps
for i=1:7
    figure(100+i), clf
    BS = BMRAresults(i).BS;
    PS = BMRAresults(i).PS;
    %plot inferred network
    df_pcolor(BS.*(PS>0.0))
    hold on
    %plot initial literature network
    [ix, iy] = find(myInteractionMatrix'==1);
    scatter(ix+0.5, iy+0.5, 'xk')
    %annotate interaction stats, use threshhold on BS
    GS = abs(BS)>0.1;
    ingoing = sum(GS');
    outgoing = sum(GS);
    text([1:nNodes]+0.5, repmat(0, nNodes,1), num2str(outgoing(:)), 'HorizontalAlignment', 'center')
    text(repmat(nNodes+1.5, nNodes,1), [1:nNodes]+0.5, num2str(ingoing(:)), 'VerticalAlignment' , 'middle' )
    %make look nice
    colormap(redbluecmap)
    set(gca, 'xTick', [1:nNodes]+0, 'XTickLabel', myNetwork, 'XTickLabelRotation', 60)
    set(gca, 'yTick', [1:nNodes]+0.5, 'yTickLabel', myNetwork)
    set(gca, 'YDir', 'reverse')
    set(gca, 'CLim', [-1 1])
    th = title(myCellLines(i));
    th.Position(1) = - 3.5;
    th.Position(2) = 0.5;
    th.HorizontalAlignment = 'left';
    th.VerticalAlignment = 'bottom';
%     %
%     filename = sprintf('../figs/BMRI_NW_%s_3', myCellLines(i));
%     print(gcf, [filename '.png'], '-dpng', '-r600')
%     %ANNOTATION:
%     fileID = fopen([filename '.txt'],'w');
%     fprintf(fileID, 'Source: /code/BMRA/bmra_analysis_heatmaps.m\n');
%     fprintf(fileID, 'Date: %s\n', datestr(now));
%     fprintf(fileID, 'Datafile: %s\n', datafile);
%     fclose(fileID);
end


%% calculate consensus network sens vs resist

displaythreshold = 0;

sortresist = 0; % if 1 then consitent pos int in resitant cells are found
                %      => figure 2
                % if 0 then consitent pos int in sensitive cells are found
                %      => figure 1

myIntMat_all = [reshape(BMRAresults(1).BS',1,[])
    reshape(BMRAresults(2).BS', 1, [])
    reshape(BMRAresults(3).BS',1,[])
    reshape(BMRAresults(4).BS',1,[])
    reshape(BMRAresults(5).BS',1,[])
    reshape(BMRAresults(6).BS',1,[])
    reshape(BMRAresults(7).BS',1,[])
    ];
myIntMat_sens = myIntMat_all(1:3,:);
myIntMat_resi = myIntMat_all(5:7,:);

%%

%Sensitive:
posInt = find(all(myIntMat_sens>0));
%posInt = find(sum(myIntMat_sens>0)>1);
posBS = mean(myIntMat_sens(:,posInt),1);
negInt = find(all(myIntMat_sens<0));
%negInt = find(sum(myIntMat_sens<0)>1);
negBS = mean(myIntMat_sens(:,negInt),1);

myConsNW = zeros(nNodes);
myConsNW(posInt) = posBS;
myConsNW(negInt) = negBS;
myConsNW = myConsNW'%.*(abs(myConsNW')>displaythreshold);
myConsNW_sens = myConsNW;

myConsNW(idx2,idx1)


%% Heatmap for consensus NW in sensitive cells
figure(1), clf
   df_pcolor(myConsNW)
    hold on
    %plot initial literature network
    [ix, iy] = find(myInteractionMatrix'==1);
    scatter(ix+0.5, iy+0.5, 'xk')
    %annotate interaction stats, use threshhold on BS
    GS = abs(myConsNW)>0.1;
    ingoing = sum(GS');
    outgoing = sum(GS);
    text([1:nNodes]+0.5, repmat(0, nNodes,1), num2str(outgoing(:)), 'HorizontalAlignment', 'center')
    text(repmat(nNodes+1.5, nNodes,1), [1:nNodes]+0.5, num2str(ingoing(:)), 'VerticalAlignment' , 'middle' )
    %make look nice
    colormap(redbluecmap(51))
    set(gca, 'xTick', [1:nNodes]+0, 'XTickLabel', myNetwork, 'XTickLabelRotation', 60)
    set(gca, 'yTick', [1:nNodes]+0.5, 'yTickLabel', myNetwork)
    set(gca, 'YDir', 'reverse')
    set(gca, 'CLim', 1*[-1 1])
    th = title('Sensitive');
    th.Position(1) = - 3.5;
    th.Position(2) = 0.5;
    th.HorizontalAlignment = 'left';
    th.VerticalAlignment = 'bottom';

%    filename = sprintf('../figs/BMRI_consNW_Sensitive', myCellLines(i));
%    print(gcf, [filename '.png'], '-dpng', '-r600')
%    %ANNOTATION:
%    fileID = fopen([filename '.txt'],'w');
%    fprintf(fileID, 'Source: /code/BMRA/bmra_analysis_heatmaps.m\n');
%    fprintf(fileID, 'Date: %s\n', datestr(now));
%    fprintf(fileID, 'Datafile: %s\n', datafile);
%    fclose(fileID);

 [i1 i2 ] = find(abs(myConsNW )~=0);
 [myNetwork(i2) myNetwork(i1)]

%% Heatmap for consensus NW in resistant cells
posInt = find(all(myIntMat_resi>0));
%posInt = find(sum(myIntMat_resi>0)>1);
posBS = mean(myIntMat_resi(:,posInt),1);
negInt = find(all(myIntMat_resi<0));
%negInt = find(sum(myIntMat_resi<0)>1);
negBS = mean(myIntMat_resi(:,negInt),1);

myConsNW = zeros(nNodes);
myConsNW(posInt) = posBS;
myConsNW(negInt) = negBS;
myConsNW = myConsNW'.*(abs(myConsNW')>displaythreshold);
myCOnsNW_resi = myConsNW;

figure(2), clf
   df_pcolor(myConsNW)
    hold on
    %plot initial literature network
    scatter(ix+0.5, iy+0.5, 'xk')
    %annotate interaction stats, use threshhold on BS
    GS = abs(myConsNW)>0.1;
    ingoing = sum(GS');
    outgoing = sum(GS);
    text([1:nNodes]+0.5, repmat(0, nNodes,1), num2str(outgoing(:)), 'HorizontalAlignment', 'center')
    text(repmat(nNodes+1.5, nNodes,1), [1:nNodes]+0.5, num2str(ingoing(:)), 'VerticalAlignment' , 'middle' )
    %make look nice
    colormap(redbluecmap)
    set(gca, 'xTick', [1:nNodes]+0, 'XTickLabel', myNetwork, 'XTickLabelRotation', 60)
    set(gca, 'yTick', [1:nNodes]+0.5, 'yTickLabel', myNetwork)
    set(gca, 'YDir', 'reverse')
    set(gca, 'CLim', [-1 1])
    th = title('Resistant');
    th.Position(1) = - 3.5;
    th.Position(2) = 0.5;
    th.HorizontalAlignment = 'left';
    th.VerticalAlignment = 'bottom';
    
%    save 'mat_bmra_analysis_ConsNW.mat' myConsNW_sens myCOnsNW_resi myInteractionMatrix myNetwork

%    filename = sprintf('../figures/BMRI_consNW_Resistant', myCellLines(i))
%    print(gcf, [filename '.png'], '-dpng', '-r600')
%    %ANNOTATION:
%    fileID = fopen([filename '.txt'],'w');
%    fprintf(fileID, 'Source: /code/BMRA/bmra_analysis_heatmaps.m\n');
%    fprintf(fileID, 'Date: %s\n', datestr(now));
%    fprintf(fileID, 'Datafile: %s\n', datafile);
%    fclose(fileID);
    
[i1 i2 ] = find(abs(myConsNW )~=0);
 [myNetwork(i2) myNetwork(i1)]
 