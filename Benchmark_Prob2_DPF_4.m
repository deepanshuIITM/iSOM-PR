%% INTERPRETABLE SELF ORGANIZING MAPS (iSOM) ASSISTED INTERACTIVE MULTI-CRITERIA DECISION-MAKING FOLLOWING PARETO RACE
%% DPF-4 PROBLEM: Benchmark problem 2 in the manuscript
% Data set for Pareto front generation and iteration-wise DM's preferred solutions are
% provided for replication purposes. It is to be noted the iSOM plots may
% change with changing the data and DM's preferred solutions. 
% NOTE: THE UNIFIED FRAMEWORK FOR ISOM-AIDED PARETO RACE IS PROVIDED
% ELSEWHERE. THIS CODE IS FOR REPLICATION PURPOSE ONLY.

%% NOTE: FOR DETAILED EXPLANATION, AT MULTIPLE PLACES PAUSE COMMAND IS GIVEN. SEE THE COMMAND WINDOW FOR APPROPRIATE ACTION (STRIKE ANY KEY TO PROCEED)

%%
clc
clear all
close all

%% UPLOAD FINITE SIZE REPRESENTATION OF NEAR PARETO FRONT
Tb = readtable('Benchmark_Prob2_DPF_4.csv');  d = table2array(Tb);  data = d(1:end,2:end);
% Assigning Objective Functions
F(:,1) = data(:,1); 
F(:,2) = data(:,2); 
F(:,3) = data(:,3); 
F(:,4) = data(:,4);

%% AVERAGE CV CALCULATION (ALL CONSTRAINT ARE NORMALIZED IN THE PROBLEM)
avg_CV = data(:,end);                               % Taking Average of all normalized constraint
F(:,5) = (avg_CV-min(avg_CV))/range(avg_CV) - 1;    % Creating an array for average CV

CV = avg_CV(find(avg_CV>max(avg_CV)-0.025*range(avg_CV)));  % Here we considered top 2.5% constraint, DM is free to change it
F_CV = F(find(avg_CV>max(avg_CV)-0.025*range(avg_CV)),:);

%% TRADE-OFF CALCULATION
sData_t = som_data_struct(F(:,1:4),'my-data','comp_names',{'f1','f2','f3','f4'});
sData_t = som_normalize(sData_t,'range');
f = sData_t.data;
% determine size of the data-set (n: #pts, m: #obj)
[n,m] = size(f);
% Set all points as original points
typedata = repmat("Original",n,1);
labels = {'f1','f2','f3','f4'};
% set neighborhood size 
neighbor = 10;
% Number of points to be reported
numRep = 3;
infty = 1e6;
for i=1:n
    % compute pair-wise distances, except distance from itself
    dist = [];
    for j=1:n
        if (i~=j)
            dist = [dist, norm(f(i,:)-f(j,:),2)];
        else dist = [dist, infty];
        end
    end
    % then, sort the distances in ascending order
    [val,id] = sort(dist,'ascend');
    % note the nearest 10 neighbor IDs
    sortid = id(1:neighbor);
    % compute tradeoff
    tradeof = zeros(1,m); % initialize trade-off to zero
    for l=1:neighbor
        count1 = 0; % #loss objectives
        count2 = 0; % #gain objectives
        sumloss = 0.0; sumgain = 0.0;
        for k=1:m
            if f(i,k) < f(sortid(l),k)
                count1 = count1 + 1;
                sumloss = sumloss + (f(sortid(l),k)-f(i,k)); % Calculation of average loss
            end
            if f(i,k) > f(sortid(l),k)
                count2 = count2 + 1;
                sumgain = sumgain + (f(i,k) - f(sortid(l),k));  % Calculation of average gain
            end
        end
        tradeoffNeigh(l) = (sumloss/count1)/(sumgain/count2);   % Trade-off Calculation
    end
    % compute the max trade-off of all neighbors
    tradeoff(i) = max(tradeoffNeigh);
end
F(:,6) = tradeoff';  % Creating an array for Trade-off value

%% Data Preperation for iSOM Initialization 
f1=max(F(:,1))-min(F(:,1)); f2=max(F(:,2))-min(F(:,2));
f3=max(F(:,3))-min(F(:,3)); f4=max(F(:,4))-min(F(:,4));

sum = f1+f2+f3+f4;
F_final = (f1/sum)*F(:,1)+(f2/sum)*F(:,2)+(f3/sum)*F(:,3)+(f4/sum)*F(:,4); % (Weighted Objective Vector for iSOM initialization, 
% see supplementary documents for additional details)

%% Making Data for iSOM Initialization
dataf = data(:,5:6); test =[dataf F_final];
test1=[dataf F(:,1)]; % Objective function 1
test2=[dataf F(:,2)]; % Objective function 2
test3=[dataf F(:,3)]; % Objective function 3
test4=[dataf F(:,4)]; % Objective function 4
test5=[dataf F(:,5)]; % Near Constraint Violation Metric (G)
test6=[dataf F(:,6)]; % Trade-off Metric (T)

test_F = [dataf F(:,1) F(:,2) F(:,3) F(:,4) F(:,5) F(:,6)]; % Dataset of Design Variables, Objectives, G, T

%%
sData_F = som_data_struct(test_F,'comp_names',{'x1','x2', 'F1','F2','F3','F4','G','T'}); sData_F = som_normalize(sData_F,'range'); % DATASET
sData  = som_data_struct(test,'comp_names',{'x1','x2','F'});   sData  = som_normalize(sData,'range');
sData1 = som_data_struct(test1,'comp_names',{'x1','x2','F1'}); sData1 = som_normalize(sData1,'range'); % Objective function 1
sData2 = som_data_struct(test2,'comp_names',{'x1','x2','F2'}); sData2 = som_normalize(sData2,'range'); % Objective function 2
sData3 = som_data_struct(test3,'comp_names',{'x1','x2','F3'}); sData3 = som_normalize(sData3,'range'); % Objective function 3
sData4 = som_data_struct(test4,'comp_names',{'x1','x2','F4'}); sData4 = som_normalize(sData4,'range'); % Objective function 4
sData5 = som_data_struct(test5,'comp_names',{'x1','x2','G'});  sData5 = som_normalize(sData5,'range'); % Near Constraint Violation Metric (G)
sData6 = som_data_struct(test6,'comp_names',{'x1','x2','T'});  sData6 = som_normalize(sData6,'range'); % Trade-off Metric (T)

%% INITIALIZING ISOM CODEBOOK VECTOR (Linear Initialization)
[sMap_F]= som_lininit(sData_F,'lattice','hexa','msize',[30,30]);
[sMap]= modifiedsom_lininit1(sData1,'lattice','hexa','msize',sMap_F.topol.msize);
% sMap.codebook(:,10) = sMap.codebook(:,10)*0; % optional, it does not effect the results
[sMap1]= sMap;[sMap2]= sMap; [sMap3]= sMap;
[sMap4]= sMap;[sMap5]= sMap; [sMap6]= sMap;

%% TRAINING ISOM FOR EACH OBJECTIVES, G, AND T
[sMap_F,sTrainF] = som_batchtrain(sMap_F,sData_F,'trainlen',200);
[sMap1,sTrain1] = modifiedsom_batchtrain(sMap1,sData1,'trainlen',200);
[sMap2,sTrain2] = modifiedsom_batchtrain(sMap2,sData2,'trainlen',200);
[sMap3,sTrain3] = modifiedsom_batchtrain(sMap3,sData3,'trainlen',200);
[sMap4,sTrain4] = modifiedsom_batchtrain(sMap4,sData4,'trainlen',200);
[sMap5,sTrain5] = modifiedsom_batchtrain(sMap5,sData5,'trainlen',200);
[sMap6,sTrain6] = modifiedsom_batchtrain(sMap6,sData6,'trainlen',200);

%% DENORMALIZING THE DATA
sMap_F = som_denormalize(sMap_F,sData_F); sData_F = som_denormalize(sData_F,'remove');
sMap1=som_denormalize(sMap1,sData1); sData1=som_denormalize(sData1,'remove');
sMap2=som_denormalize(sMap2,sData2); sData2=som_denormalize(sData2,'remove');
sMap3 = som_denormalize(sMap3,sData3); sData3 = som_denormalize(sData3,'remove');
sMap4 = som_denormalize(sMap4,sData4); sData4 = som_denormalize(sData4,'remove');
sMap5 = som_denormalize(sMap5,sData5); sData5 = som_denormalize(sData5,'remove');
sMap6 = som_denormalize(sMap6,sData6); sData6 = som_denormalize(sData6,'remove');

%% ARRANGING THE INDIVIDUALLY TRAINED OBJECTIVES AND PREPARING FOR PLOTTING U-MATRIX
sMap_umatrix = sMap_F;
sMap_umatrix.codebook(:,1:3) = sMap1.codebook(:,1:3);
sMap_umatrix.codebook(:,4) = sMap2.codebook(:,3);
sMap_umatrix.codebook(:,5) = sMap3.codebook(:,3);
sMap_umatrix.codebook(:,6) = sMap4.codebook(:,3);
sMap_umatrix.codebook(:,7) = sMap5.codebook(:,3);
sMap_umatrix.codebook(:,8) = sMap6.codebook(:,3);

%% iSOM PLOTS (U Matrix and Component Planes )
figure(1) 
som_show(sMap_umatrix,'umat',3:6,'comp','all');

%% TRAINED iSOM GRIDS IN OBJECTIVE FUNCTION SPACE   
figure(2)
som_grid(sMap_umatrix,'coord',sMap_umatrix.codebook(:,3:5),'label',sMap_umatrix.labels,'labelcolor','c','labelsize',10, 'marker','o','MarkerColor','k'...
    ,'MarkerSize',7,'linecolor', 'k');
hold on,
scatter3(F(:,1),F(:,2),F(:,3),20,'ro','filled');
title('Trained iSOM grid in objective Space')
xlabel('F1')
ylabel('F2')
zlabel('F3')

%% MAP SIZE / TOPOLOGY 
map_s = sMap_F.topol.msize;

%% IDENTIFYING PARETO DISJOINT BOUNDARY
dmat = som_normalize(som_dmat(sMap_umatrix.codebook(:,3:6)),'range');
h_d = zeros(map_s(1)*map_s(2),1);
for i = 1:size(h_d ,1)
    if dmat(i)> 0.485  % This value is calculated iteratively based on distance of iSOM nodes, depends on  Pareto Front data
        h_d(i,1) = 1;
    end
end

%% HIGHLIGHT PARETO DISJOINT BOUNDARY
figure(1) 
som_show(sMap_umatrix,'umat',3:6,'comp','all');
som_show_add('hit',h_d,'Markersize',1,'MarkerColor','none','EdgeColor','k','Subplot',2:8);
sgtitle('Highlighting Pareto Disjoint Boundary using Black edge color iSOM nodes'); pause(5)

%% BMU SELECTION AND VISUALIZATION OF AVG_CV
BMUsCV = som_bmus41(sMap_umatrix.codebook(:,3:end-1), F_CV);
hcv = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUsCV,1)
    n = BMUsCV(i);
    hcv(n,1)= 0.5+abs(CV(i));
end

%% HIGHLIGHT PARETO DISJOINT BOUNDARY
figure(1) 
som_show(sMap_umatrix,'umat',3:6,'comp','all');
som_show_add('hit',h_d,'Markersize',1,'MarkerColor','none','EdgeColor','k','Subplot',2:8);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:8);
sgtitle('Highlighting near CV points using Gray face color iSOM nodes')

%% Pareto Solutions
disp('Strike any key to supply Pareto Optimal Solution corresponding to start point...')
pause
start_p = [75.000 75.000 75.000 75.000];  % START POINT USED IN MANUSCRIPT (DM can provide first PO solution here)

%% BMU SELECTION OF START POINT
BMUs = som_bmus41(sMap_umatrix.codebook(:,3:end-2), start_p);
h = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs,1)
    n = BMUs(i);
    h(n,1)= 2;
end

%% HIGHLIGHTING START POINT
figure(1) 
som_show(sMap_umatrix,'umat',3:6,'comp','all');
som_show_add('hit',h_d,'Markersize',1,'MarkerColor','none','EdgeColor','k','Subplot',2:8);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:8);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:9);
sgtitle('Highlighting start point using black face color iSOM nodes')

%% DM PROVIDES THE FIRST PREFERRED SOLUTIONS OBTAINED FROM PARETO RACE
disp('Strike any key to plot first preferred solution...')
pause
iter1 = [[90.671 72.181 74.763 50.778]; [76.932  89.602 52.048 69.420]; [51.150 74.132 72.604  90.281]; [77.771 52.276  89.493 68.508]];


%% BMU SELECTION OF FIRST PREFERRED SOLUTION
BMUs1 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1);
h1 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs1,1)
    n = BMUs1(i);
    h1(n,1)= 2;
end

%% HIGHLIGHTING FIRST PREFERRED SOLUTION
figure(1) 
som_show(sMap_umatrix,'umat',3:6,'comp','all');
som_show_add('hit',h_d,'Markersize',1,'MarkerColor','none','EdgeColor','k','Subplot',2:8);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:8);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:9);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:9);
sgtitle('Highlighting first preferred solutions using green face color iSOM node')

%% DM PROVIDES THE SECOND PREFERRED SOLUTIONS OBTAINED FROM PARETO RACE
disp('Strike any key to plot second preferred solution...')
pause
iter2 = [[132.235 95.768 91.689 9.602]; [95.984 129.952 12.549 88.498]; [10.317 92.570 93.434 131.121]; [96.128 11.066 131.189 89.958]];

%% BMU SELECTION OF SECOND PREFERRED SOLUTION
BMUs2 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter2);
h2 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs2,1)
    n = BMUs2(i);
    h2(n,1)= 2;
end

%% HIGHLIGHTING SECOND PREFERRED SOLUTION
figure(1) 
som_show(sMap_umatrix,'umat',3:6,'comp','all');
som_show_add('hit',h_d,'Markersize',1,'MarkerColor','none','EdgeColor','k','Subplot',2:8);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:8);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:9);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:9);
som_show_add('hit',h2,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:9);
sgtitle('Highlighting second preferred solutions using green face color iSOM node')

%% HIGHLIGHTING TURN 
disp('Strike any key to highlight Turn...')
pause
iter3 = [[132.235 95.768 91.689 9.602]; [95.984 129.952 12.549 88.498]];

%% BMU SELECTION OF TURNING
BMUs3 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter3);
h3 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs3,1)
    n = BMUs3(i);
    h3(n,1)= 2;
end
%% HIGHLIGHTING THE SECOND PREFERRED SOLUTIONs FOR TURN
figure(1) 
som_show(sMap_umatrix,'umat',3:6,'comp','all');
som_show_add('hit',h_d,'Markersize',1,'MarkerColor','none','EdgeColor','k','Subplot',2:8);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:8);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:9);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:9);
som_show_add('hit',h2,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:9);
som_show_add('hit',h3,'Markersize',1,'MarkerColor','m','EdgeColor','m','Subplot',2:9);
sgtitle('Highlighting second preferred solutions for turn using magenta face color iSOM node')

%% HIGHLIGHTING THE THIRD AND FINAL SOLUTION
disp('Strike any key to highlight the Third and Final Solution...')
pause
iter4= [[130.116 101.393 83.248 16.752]; [101.643 129.069 18.362 81.638]; [0.0000 100.00 100.00 141.421];[100.00 0.0000 141.421 100.00]]; 

%% BMU SELECTION OF TURNING
BMUs4 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter4);
h4 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs4,1)
    n = BMUs4(i);
    h4(n,1)= 2;
end
%% HIGHLIGHTING THE THIRD PREFERRED SOLUTIONs FOR TURN
figure(1) 
som_show(sMap_umatrix,'umat',3:6,'comp','all');
som_show_add('hit',h_d,'Markersize',1,'MarkerColor','none','EdgeColor','k','Subplot',2:8);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:8);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:9);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:9);
som_show_add('hit',h2,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:9);
som_show_add('hit',h3,'Markersize',1,'MarkerColor','m','EdgeColor','m','Subplot',2:9);
som_show_add('hit',h4,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:9);
som_show_add('hit',h4,'Markersize',1,'MarkerColor','r','EdgeColor','r','Subplot',2:9);
sgtitle('Highlighting second preferred solutions for turn using red face color iSOM node')

%% Scatter Plot of Objectives ( G and T metric is not plotted here)
iterx= [[90.671 72.181 74.763 50.778] 
        [132.235 95.768 91.689 9.602]
        [130.116 101.393 83.248 16.752]        

        [76.932  89.602 52.048 69.420]
        [95.984 129.952 12.549 88.498]
        [101.643 129.069 18.362 81.638]

        [51.150 74.132 72.604  90.281]
        [10.317 92.570 93.434 131.121]
        [0.0000 100.00 100.00 141.421]

        [77.771 52.276  89.493 68.508]
        [96.128 11.066 131.189 89.958]
        [100.00 0.0000 141.421 100.00]
         ]; 

figure(3)
scatter3(data(:,1), data(:,2), data(:,3), 2, data(:,4), 'MarkerFaceColor', 'Flat','Marker','o'); 
hold on;
colormap('parula')
scatter3(start_p(:,1),start_p(:,2),start_p(:,3),50,'dk','filled');
scatter3(iterx(:,1), iterx(:,2), iterx(:,3),30, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'k' ,'Marker','s');
scatter3(iterx([2,5],1), iterx([2,5],2), iterx([2,5],3),32, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'k' ,'Marker','s');
scatter3(iterx([3 6 9 12],1), iterx([3 6 9 12],2), iterx([3 6 9 12],3),30, 'r','filled','Marker','s');

colorbar;
caxis([min(data(:,4)), max(data(:,4))])

legend('Pareto Front','Start Point','Constant Speed','Turn','Final Solution','Location','best')
xlabel('$f_1$','Interpreter','latex','FontSize',15)
ylabel('$f_2$','Interpreter','latex','FontSize',15)
zlabel('$f_3$','Interpreter','latex','FontSize',15)
sgtitle('PR Solutions in Objective Space (G and T plots are not shown)')