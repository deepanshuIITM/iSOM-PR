%%
clc;
close all;
clear all;

%%
Tb = readtable('Article_New.csv');  
d = table2array(Tb);
data = d(:,2:end);
F(:,1) = data(:,1); F(:,2) = data(:,2); F(:,3) = data(:,3); 

start_p1 = [-0.164  -3.083  -14.695]; 
iter1 = [[-0.067  -2.896 -16.948]
         [0.314  -2.256 -22.516]
         [1.139  -0.979 -30.697]           
         [1.866   0.116 -36.814]
         [2.443   0.979 -41.521]
         [3.602   1.144 -47.069]
        ];

%% Plot the solution points.
figure(1)
sz = 5;
scatter3(data(:,1), data(:,2), data(:,3), 2, 'b', 'MarkerFaceColor', 'Flat','Marker','o'); hold on;
scatter3(start_p1(:,1), start_p1(:,2), start_p1(:,3), 30, 'k', 'MarkerFaceColor', 'Flat','Marker','d');
scatter3(iter1(:,1), iter1(:,2), iter1(:,3), 30, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'k' ,'Marker','s');
scatter3(iter1(end,1), iter1(end,2), iter1(end,3), 30, 'r', 'MarkerFaceColor', 'Flat','Marker','s');
xlabel('$f_1$','Interpreter','latex','FontSize',15)
ylabel('$f_2$','Interpreter','latex','FontSize',15)
zlabel('$f_3$','Interpreter','latex','FontSize',15)
legend('Pareto Front','Start Point', 'Constant Speed', 'Final Sol')

%%
f1=max(F(:,1))-min(F(:,1)); f2=max(F(:,2))-min(F(:,2)); f3=max(F(:,3))-min(F(:,3));  
sum = f1+f2+f3;
F_final = (f1/sum)*F(:,1)+(f2/sum)*F(:,2)+(f3/sum)*F(:,3);

%% average CV calculation 
avg_CV = mean(data(:,4:6)');
F(:,4) = avg_CV;
CV = avg_CV(find(avg_CV>max(avg_CV)-0.05*range(avg_CV)));
F_CV = F(find(avg_CV>max(avg_CV)-0.05*range(avg_CV)),:);
[val,idx] = max(data(:,4:6)');
idx1 = idx(find(avg_CV>max(avg_CV)-0.05*range(avg_CV)));

%% Trade-off Calculation
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
                sumloss = sumloss + (f(sortid(l),k)-f(i,k));
            end
            if f(i,k) > f(sortid(l),k)
                count2 = count2 + 1;
                sumgain = sumgain + (f(i,k) - f(sortid(l),k));
            end
        end
        tradeoffNeigh(l) = (sumloss/count1)/(sumgain/count2);
    end
    % compute the max trade-off of all neighbors
    tradeoff(i) = max(tradeoffNeigh);
end
F(:,5) = tradeoff';

%% Polynomial 
%csvwrite('Clutch_plot.csv', F)

% T = array2table(F);
% T.Properties.VariableNames(1:5) = {'F1','F2','F3','G','T'};
% writetable(T,'Article_plot.csv')

%% Making DataStruct for SOM Training
test =[data(:,7:end) F_final];
test1=[data(:,7:end) F(:,1)]; test2=[data(:,7:end) F(:,2)];
test3=[data(:,7:end) F(:,3)]; test4=[data(:,7:end) F(:,4)]; test5=[data(:,7:end) F(:,5)];
test_F = [data(:,7:end) F(:,1) F(:,2) F(:,3) F(:,4) F(:,5)];

%%
sData_F = som_data_struct(test_F,'comp_names',{'x1','x2','F1','F2','F3','G','T'});
sData_F = som_normalize(sData_F,'range');

sData = som_data_struct(test,'comp_names',{'x1','x2','F'});    sData = som_normalize(sData,'range');
sData1 = som_data_struct(test1,'comp_names',{'x1','x2','F1'}); sData1 = som_normalize(sData1,'range');
sData2 = som_data_struct(test2,'comp_names',{'x1','x2','F2'}); sData2 = som_normalize(sData2,'range');
sData3 = som_data_struct(test3,'comp_names',{'x1','x2','F3'}); sData3 = som_normalize(sData3,'range');
sData4 = som_data_struct(test4,'comp_names',{'x1','x2','G'}); sData4 = som_normalize(sData4,'range');
sData5 = som_data_struct(test5,'comp_names',{'x1','x2','T'}); sData5 = som_normalize(sData5,'range');

%% Initializing SOM Map Codebook Vectors (Linear Initialization)
[sMap_F]= som_lininit(sData_F,'lattice','hexa','msize',[20,20]);
[sMap]= modifiedsom_lininit(sData,'lattice','hexa','msize',[20,20]);
% sMap.codebook(:,6) = sMap.codebook(:,6)*0; % optional, it does not effect the results
[sMap1]= sMap;[sMap2]= sMap;[sMap3]= sMap;[sMap4]= sMap;[sMap5]= sMap;[sMap6]= sMap;

%% Training SOM
[sMap_F,sTrainF] = som_batchtrain(sMap_F,sData_F,'sample_order','ordered','trainlen',500,'rad_init',1.25, 'rad_fin', 0.5);
[sMap1,sTrain1] = modifiedsom_batchtrain(sMap1,sData1,'sample_order','ordered','radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500);
[sMap2,sTrain2] = modifiedsom_batchtrain(sMap2,sData2,'sample_order','ordered','radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500);
[sMap3,sTrain3] = modifiedsom_batchtrain(sMap3,sData3,'sample_order','ordered','radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500);
[sMap4,sTrain4] = modifiedsom_batchtrain(sMap4,sData4,'sample_order','ordered','radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500);
[sMap5,sTrain5] = modifiedsom_batchtrain(sMap5,sData5,'sample_order','ordered','radius_ini', 1.0, 'radius_fin', 0.75,'trainlen',500);

%% Denormalizing the data
sMap_F = som_denormalize(sMap_F,sData_F); sData_F = som_denormalize(sData_F,'remove');
sMap1=som_denormalize(sMap1,sData1);   sData1=som_denormalize(sData1,'remove');
sMap2=som_denormalize(sMap2,sData2);   sData2=som_denormalize(sData2,'remove');
sMap3 = som_denormalize(sMap3,sData3); sData3 = som_denormalize(sData3,'remove');
sMap4 = som_denormalize(sMap4,sData4); sData4 = som_denormalize(sData4,'remove');
sMap5 = som_denormalize(sMap5,sData5); sData5 = som_denormalize(sData5,'remove');

%%
sMap_umatrix = sMap_F;
sMap_umatrix.codebook(:,1:3) = sMap1.codebook(:,1:3);
sMap_umatrix.codebook(:,4) = sMap2.codebook(:,3);
sMap_umatrix.codebook(:,5) = sMap3.codebook(:,3);
sMap_umatrix.codebook(:,6) = sMap4.codebook(:,3);
sMap_umatrix.codebook(:,7) = sMap5.codebook(:,3);

%% Visualization of SOM results( U Matrix and Component Planes )
figure(2) 
som_show(sMap_F);

figure(3) 
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);

%% iSOM Grid in function space  
figure(4)
som_grid(sMap_umatrix,'coord',sMap_umatrix.codebook(:,[3 4 5]),'label',sMap_umatrix.labels,'labelcolor','c','labelsize',10, 'marker','o','MarkerColor','k'...
    ,'MarkerSize',7,'linecolor', 'k');
hold on, scatter3(F(:,1),F(:,2),F(:,3),20,'ro','filled');
xlabel('F1')
ylabel('F2')
zlabel('F3')

%% cSOM Grid in function space 
figure(5)
som_grid(sMap_F,'coord',sMap_F.codebook(:,[3 4 5]),'label',sMap_F.labels,'labelcolor','c','labelsize',10, 'marker','o','MarkerColor','k'...
    ,'MarkerSize',7,'linecolor', 'k');
hold on, scatter3(F(:,1),F(:,2),F(:,3),20,'ro','filled');
xlabel('F1')
ylabel('F2')
zlabel('F3')

%% BMU selection and visualization of starting point
map_s = [20 20];
BMUs1 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), start_p1);
h1 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs1,1)
    n = BMUs1(i);
    h1(n,1)= 2;
end

%% BMU selection and visualization of avg_CV
BMUs2 = som_bmus41(sMap_umatrix.codebook(:,3:end-1), F_CV);
h2 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs2,1)
    n = BMUs2(i);
    h2(n,1)= 0.5+abs(CV(i));
end

%% BMU selection and visualization of sol 1 
BMUs31 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1(1,:));
h31 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs31(:,1),1)
    n = BMUs31(i);
    h31(n,1)= h31(n,1)+1;
end

%% BMU selection and visualization of sol 2 
BMUs32 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1(2,:));
h32 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs32(:,1),1)
    n = BMUs32(i);
    h32(n,1)= h32(n,1)+1;
end

%% BMU selection and visualization of sol 3
BMUs33 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1(3,:));
h33 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs33(:,1),1)
    n = BMUs33(i);
    h33(n,1)= h33(n,1)+1;
end
%% BMU selection and visualization of sol 4
BMUs34 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1(4,:));
h34 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs34(:,1),1)
    n = BMUs34(i);
    h34(n,1)= h34(n,1)+1;
end

%% BMU selection and visualization of sol 5
BMUs35 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1(5,:));
h35 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs35(:,1),1)
    n = BMUs35(i);
    h35(n,1)= h35(n,1)+1;
end

%% BMU selection and visualization of sol 6 
BMUs36 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1(6,:));
h36 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs36(:,1),1)
    n = BMUs36(i);
    h36(n,1)= h36(n,1)+1;
end

%% BMU selection for final solution
BMUs4 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1(end,:));
h4 = zeros(map_s(1)*map_s(2),1);
j = 1;
for i = 1:size(BMUs4(:,1),1)
    n = BMUs4(i);
    h4(n,1)= h4(n,1)+1;
end

%% BMU selection and visualization of iter 2 
BMUs5 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1([4 6],:));
h5 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs5(:,1),1)
    n = BMUs5(i);
    h5(n,1)= h5(n,1)+1;
end 

%% BMU selection and visualization of iter 3 
BMUs6 = som_bmus41(sMap_umatrix.codebook(:,3:end-2), iter1([2 6],:));
h6 = zeros(map_s(1)*map_s(2),1);

for i = 1:size(BMUs6(:,1),1)
    n = BMUs6(i);
    h6(n,1)= h6(n,1)+1;
end 

%% Final plots
figure(8)
som_show(sMap_umatrix,'umat',3:7,'comp',3:7);
% som_show(sMap2,'umat',1:4,'comp',4,'bar','none');
saveas(figure(8),sprintf('Figure_%d.png',1))
som_show_add('hit',h2,'Markersize',1,'MarkerColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1],'Subplot',2:6);
saveas(figure(8),sprintf('Figure_%d.png',2))
som_show_add('hit',h1,'Markersize',1,'MarkerColor','k','EdgeColor','k','Subplot',2:6);
saveas(figure(8),sprintf('Figure_%d.png',3))
som_show_add('hit',h31,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
saveas(figure(8),sprintf('Figure_%d.png',4))
som_show_add('hit',h32,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
saveas(figure(8),sprintf('Figure_%d.png',5))
som_show_add('hit',h33,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
saveas(figure(8),sprintf('Figure_%d.png',6))
som_show_add('hit',h34,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
saveas(figure(8),sprintf('Figure_%d.png',7))
som_show_add('hit',h35,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
saveas(figure(8),sprintf('Figure_%d.png',8))
som_show_add('hit',h36,'Markersize',1,'MarkerColor','g','EdgeColor','g','Subplot',2:6);
saveas(figure(8),sprintf('Figure_%d.png',9))
som_show_add('hit',h4,'Markersize',1,'MarkerColor','r','EdgeColor','r','Subplot',2:6);
saveas(figure(8),sprintf('Figure_%d.png',10))
