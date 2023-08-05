%% INTERPRETABLE SELF ORGANIZING MAPS (iSOM) ASSISTED INTERACTIVE MULTI-CRITERIA DECISION-MAKING FOLLOWING PARETO RACE
% ARTICLE PROBLEM: Benchmark problem 1 in the manuscript
% Data set for Pareto front generation and iteration-wise DM's preferred solutions are
% provided for replication purposes. It is to be noted the iSOM plots may
% change with changing the data and DM's preferred solutions. 
% NOTE: THE UNIFIED FRAMEWORK FOR ISOM-AIDED PARETO RACE IS PROVIDED
% ELSEWHERE. THIS CODE IS FOR REPLICATION PURPOSE ONLY.

%% NOTE: FOR DETAILED EXPLANATION, AT MULTIPLE PLACES PAUSE COMMAND IS GIVEN. SEE THE COMMAND WINDOW FOR APPROPRIATE ACTION (STRIKE ANY KEY TO PROCEED)

%% 
clc;
close all;
clear all;

%% UPLOAD FINITE SIZE REPRESENTATION OF NEAR PARETO FRONT
Tb = readtable('Benchmark_Prob1.csv');  
d = table2array(Tb);
data = d(:,2:end);
% Assigning Objective Functions
F(:,1) = data(:,1); F(:,2) = data(:,2); F(:,3) = data(:,3); 

%% AVERAGE CV CALCULATION (ALL CONSTRAINT ARE NORMALIZED IN THE PROBLEM)
avg_CV = mean(data(:,4:6)');  % Taking Average of all normalized constraint
F(:,4) = avg_CV;              % Creating an array for average CV
CV = avg_CV(find(avg_CV>max(avg_CV)-0.05*range(avg_CV))); % Here we considered top 5% constraint, DM is free to change it
F_CV = F(find(avg_CV>max(avg_CV)-0.05*range(avg_CV)),:);
[val,idx] = max(data(:,4:6)');
idx1 = idx(find(avg_CV>max(avg_CV)-0.05*range(avg_CV)));

%% TRADE-OFF CALCULATION
sData_t = som_data_struct(F(:,1:3),'my-data','comp_names',{'f1','f2','f3'});
sData_t = som_normalize(sData_t,'range');
f = sData_t.data;
% determine size of the data-set (n: #pts, m: #obj)
[n,m] = size(f);
% Set all points as original points
typedata = repmat("Original",n,1);
labels = {'f1','f2','f3'};
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
                sumloss = sumloss + (f(sortid(l),k)-f(i,k));    % Calculation of average loss
            end
            if f(i,k) > f(sortid(l),k)
                count2 = count2 + 1;
                sumgain = sumgain + (f(i,k) - f(sortid(l),k));  % Calculation of average gain
            end
        end
        tradeoffNeigh(l) = (sumloss/count1)/(sumgain/count2);  % Trade-off Calculation 
    end
    % compute the max trade-off of all neighbors
    tradeoff(i) = max(tradeoffNeigh);
end
F(:,5) = tradeoff';  % Creating an array for Trade-off value

%% Data Preperation for iSOM Initialization 
f1=max(F(:,1))-min(F(:,1)); f2=max(F(:,2))-min(F(:,2)); f3=max(F(:,3))-min(F(:,3));  
sum = f1+f2+f3;
F_final = (f1/sum)*F(:,1)+(f2/sum)*F(:,2)+(f3/sum)*F(:,3);  % (Weighted Objective Vector for iSOM initialization, 
% see supplementary documents for additional details)

%% Making Data for iSOM Initialization
test =[data(:,7:end) F_final];
test1=[data(:,7:end) F(:,1)]; % Objective function 1
test2=[data(:,7:end) F(:,2)]; % Objective function 2
test3=[data(:,7:end) F(:,3)]; % Objective function 3
test4=[data(:,7:end) F(:,4)]; % Near Constraint Violation Metric (G)
test5=[data(:,7:end) F(:,5)]; % Trade-off Metric (T)
test_F = [data(:,7:end) F(:,1) F(:,2) F(:,3) F(:,4) F(:,5)]; % Dataset of Design Variables, Objectives, G, T

%% NORMALIZING DATA FOR INITIALIZATION 
sData_F = som_data_struct(test_F,'comp_names',{'x1','x2','F1','F2','F3','G','T'}); sData_F = som_normalize(sData_F,'range'); % DATASET

sData  = som_data_struct(test,'comp_names',{'x1','x2','F'});    sData = som_normalize(sData,'range');
sData1 = som_data_struct(test1,'comp_names',{'x1','x2','F1'}); sData1 = som_normalize(sData1,'range');  % Objective function 1
sData2 = som_data_struct(test2,'comp_names',{'x1','x2','F2'}); sData2 = som_normalize(sData2,'range');  % Objective function 2
sData3 = som_data_struct(test3,'comp_names',{'x1','x2','F3'}); sData3 = som_normalize(sData3,'range');  % Objective function 3
sData4 = som_data_struct(test4,'comp_names',{'x1','x2','G'});  sData4 = som_normalize(sData4,'range');  % Near Constraint Violation Metric (G)
sData5 = som_data_struct(test5,'comp_names',{'x1','x2','T'});  sData5 = som_normalize(sData5,'range');  % Trade-off Metric (T)

%% INITIALIZING ISOM CODEBOOK VECTOR (Linear Initialization)
[sMap_F]= som_lininit(sData_F,'lattice','hexa','msize',[20,20]);

[sMap]= modifiedsom_lininit(sData,'lattice','hexa','msize',[20,20]); % map size can be decided by iSOM itself, we used [20,20] here, DM can change it.
[sMap1]= sMap;[sMap2]= sMap;[sMap3]= sMap;[sMap4]= sMap;[sMap5]= sMap;[sMap6]= sMap;

%% TRAINING ISOM FOR EACH OBJECTIVES, G, AND T
[sMap_F,sTrainF] = som_batchtrain(sMap_F,sData_F,'trainlen',500,'rad_init',1.25, 'rad_fin', 0.5);
[sMap1,sTrain1] = modifiedsom_batchtrain(sMap1,sData1,'radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500); % OBJECTIVE 1
[sMap2,sTrain2] = modifiedsom_batchtrain(sMap2,sData2,'radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500); % OBJECTIVE 2
[sMap3,sTrain3] = modifiedsom_batchtrain(sMap3,sData3,'radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500); % OBJECTIVE 3
[sMap4,sTrain4] = modifiedsom_batchtrain(sMap4,sData4,'radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500); % G METRIC
[sMap5,sTrain5] = modifiedsom_batchtrain(sMap5,sData5,'radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500); % T METRIC

%% DENORMALIZING THE DATA
sMap_F = som_denormalize(sMap_F,sData_F); sData_F = som_denormalize(sData_F,'remove');
sMap1 = som_denormalize(sMap1,sData1);    sData1 = som_denormalize(sData1,'remove');
sMap2 = som_denormalize(sMap2,sData2);    sData2 = som_denormalize(sData2,'remove');
sMap3 = som_denormalize(sMap3,sData3);    sData3 = som_denormalize(sData3,'remove');
sMap4 = som_denormalize(sMap4,sData4);    sData4 = som_denormalize(sData4,'remove');
sMap5 = som_denormalize(sMap5,sData5);    sData5 = som_denormalize(sData5,'remove');

%% ARRANGING THE INDIVIDUALLY TRAINED OBJECTIVES AND PREPARING FOR PLOTTING U-MATRIX 
sMap_umatrix = sMap_F;
sMap_umatrix.codebook(:,1:3) = sMap1.codebook(:,1:3);
sMap_umatrix.codebook(:,4) = sMap2.codebook(:,3);
sMap_umatrix.codebook(:,5) = sMap3.codebook(:,3);
sMap_umatrix.codebook(:,6) = sMap4.codebook(:,3);
sMap_umatrix.codebook(:,7) = sMap5.codebook(:,3);

%% iSOM PLOTS (U Matrix and Component Planes )
figure(1) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7); 
sgtitle('Component planes of U-matrix,objective functions, G, and T')

%% TRAINED iSOM GRIDS IN OBJECTIVE FUNCTION SPACE  
figure(2)
som_grid(sMap_umatrix,'coord',sMap_umatrix.codebook(:,[3 4 5]),'label',sMap_umatrix.labels,'labelcolor','c','labelsize',10, 'marker','o','MarkerColor','k'...
    ,'MarkerSize',7,'linecolor', 'k');
hold on, scatter3(F(:,1),F(:,2),F(:,3),20,'ro','filled');
title('Trained iSOM grid in objective Space')
xlabel('F1')
ylabel('F2')
zlabel('F3')

%% BMU SELECTION AND VISUALIZATION OF AVG_CV
map_s = [20 20];
BMUsCV = som_bmus41(sMap_umatrix.codebook(:,3:end-1), F_CV);
hcv = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUsCV ,1)
    n = BMUsCV (i);
    hcv(n,1)= 0.5+abs(CV(i));
end

%% HIGHLIGHTING NEAR CONSTRAINT POINTS
figure(1) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
sgtitle('Highlighting near constraint points using Grey face color iSOM node')

%% DM PROVDES START POINT 
disp('Strike any key to supply Pareto Optimal Solution corresponding to start point...')
pause 
start_p = [-0.164  -3.083  -14.695];   % START POINT USED IN MANUSCRIPT (DM can provide first PO solution here)

%% BMU SELECTION OF START POINT
map_s = [20 20];
BMUs = som_bmus41(sMap_umatrix.codebook(:,3:end-2), start_p);
h = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs,1)
    n = BMUs(i);
    h(n,1)= 2;
end

%% HIGHLIGHTING START POINT
figure(1) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:6);
sgtitle('Highlighting Start Point using Black face color iSOM node')

%% DM PROVIDES THE FIRST PREFERRED SOLUTION OBTAINED FROM PARETO RACE
disp('Strike any key to plot first preferred solution...')
pause

iter1 = [-0.067  -2.896 -16.948];

%% BMU SELECTION OF FIRST PREFERRED SOLUTION
map_s = [20 20];
BMUs1 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1);
h1 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs1,1)
    n = BMUs1(i);
    h1(n,1)= 2;
end

%% HIGHLIGHTING FIRST PREFERRED SOLUTION
figure(1) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:6);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
sgtitle('Highlighting preferred solution using green color iSOM node')


%% DM PROVIDES THE SECOND PREFERRED SOLUTION OBTAINED FROM PARETO RACE
disp('Strike any key to plot second preferred solution...')
pause

iter2 = [0.314  -2.256 -22.516];

%% BMU SELECTION OF SECOND PREFERRED SOLUTION
map_s = [20 20];
BMUs2 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter2);
h2 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs2,1)
    n = BMUs2(i);
    h2(n,1)= 2;
end

%% HIGHLIGHTING SECOND PREFERRED SOLUTION
figure(1) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:6);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h2,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
sgtitle('Highlighting next preferred solution using green color iSOM node')

%% DM PROVIDES THE NEXT PREFERRED SOLUTION OBTAINED FROM PARETO RACE
disp('Strike any key to plot next preferred solution...')
pause

iter3 = [1.139  -0.979 -30.697];

%% BMU SELECTION OF NEXT PREFERRED SOLUTION
map_s = [20 20];
BMUs3 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter3);
h3 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs3,1)
    n = BMUs3(i);
    h3(n,1)= 2;
end

%% HIGHLIGHTING NEXT PREFERRED SOLUTION
figure(1) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:6);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h2,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h3,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
sgtitle('Highlighting next preferred solution using green color iSOM node')

%% DM PROVIDES THE NEXT PREFERRED SOLUTION OBTAINED FROM PARETO RACE
disp('Strike any key to plot next preferred solution...')
pause

iter4 = [1.866   0.116 -36.814];

%% BMU SELECTION OF NEXT PREFERRED SOLUTION
map_s = [20 20];
BMUs4 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter4);
h4 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs4,1)
    n = BMUs4(i);
    h4(n,1)= 2;
end

%% HIGHLIGHTING NEXT PREFERRED SOLUTION
figure(1) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:6);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h2,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h3,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h4,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
sgtitle('Highlighting next preferred solution using green color iSOM node')

%% DM PROVIDES THE NEXT PREFERRED SOLUTION OBTAINED FROM PARETO RACE
disp('Strike any key to plot next preferred solution...')
pause

iter5 = [2.443   0.979 -41.521];

%% BMU SELECTION OF NEXT PREFERRED SOLUTION
map_s = [20 20];
BMUs5 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter5);
h5 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs5,1)
    n = BMUs5(i);
    h5(n,1)= 2;
end

%% HIGHLIGHTING NEXT PREFERRED SOLUTION
figure(1) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:6);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h2,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h3,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h4,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h5,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
sgtitle('Highlighting next preferred solution using green color iSOM node')

%% DM PROVIDES THE NEXT PREFERRED SOLUTION OBTAINED FROM PARETO RACE
disp('Strike any key to plot next preferred solution...')
pause

iter6 = [3.602   1.144 -47.069];

%% BMU SELECTION OF NEXT PREFERRED SOLUTION
map_s = [20 20];
BMUs6 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter6);
h6 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs6,1)
    n = BMUs6(i);
    h6(n,1)= 2;
end

%% HIGHLIGHTING NEXT PREFERRED SOLUTION
figure(1) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:6);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h2,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h3,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h4,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h5,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h6,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
sgtitle('Highlighting next preferred solution using green color iSOM node')

%% PLOTTING ALL PAST SOLUTIONS ALONG WITH FINAL SOLUTION
figure(1)
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
som_show_add('hit',hcv,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
som_show_add('hit',h,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:6);
som_show_add('hit',h1,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h2,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h3,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h4,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h5,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
som_show_add('hit',h6,'Markersize',1,'MarkerColor','r','EdgeColor','r','Subplot',2:6);
sgtitle('Highlighting final solution using red color iSOM node')

%% Scatter Plot of Objectives ( G and T metric is not plotted here)
iteration = [[-0.067  -2.896 -16.948]
         [0.314  -2.256 -22.516]
         [1.139  -0.979 -30.697]           
         [1.866   0.116 -36.814]
         [2.443   0.979 -41.521]
         [3.602   1.144 -47.069]
        ];

figure(3)
scatter3(data(:,1), data(:,2), data(:,3), 2, 'b', 'MarkerFaceColor', 'Flat','Marker','o'); hold on;
scatter3(start_p(:,1), start_p(:,2), start_p(:,3), 30, 'k', 'MarkerFaceColor', 'Flat','Marker','d');
scatter3(iteration(:,1), iteration(:,2), iteration(:,3), 30, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'k' ,'Marker','s');
scatter3(iteration(end,1), iteration(end,2), iteration(end,3), 30, 'r', 'MarkerFaceColor', 'Flat','Marker','s');
xlabel('$f_1$','Interpreter','latex','FontSize',15)
ylabel('$f_2$','Interpreter','latex','FontSize',15)
zlabel('$f_3$','Interpreter','latex','FontSize',15)
legend('Pareto Front','Start Point', 'Constant Speed', 'Final Sol')
sgtitle('PR Solutions in Objective Space (G and T plots are not shown)')
