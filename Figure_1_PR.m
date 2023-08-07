%% INTERPRETABLE SELF ORGANIZING MAPS (iSOM) ASSISTED INTERACTIVE MULTI-CRITERIA DECISION-MAKING FOLLOWING PARETO RACE
%% FIGURE 1 DETAILS: START POINT, REFERENCE POINTS, REFERENCE DIRECTIONS,
% AND PREFERRED SOLUTIONS FOR FIGURE 1. THE CODE FOR GENERATING PARETO RACE
% SOLUTIONS ARE GIVEN IN .ipynb FILE

%%
clc
clear all
close all

%% IMPORTING DATA FOR PARETO FRONT 
Tb = readtable('Figure_1_PR.csv');  
d = table2array(Tb);
F = d(2:end,2:end);

%% IDEAL POINT AND NADIR POINT
ideal = [0, 0, 0];
nadir = [1, 1, 1];

%% START POINT
r_0 = [0.2, 0.2, 1.25];                                                     %USER DEFINED START POINT
asp = [1.0, 0.0, 0.0];                                                      %CLASSIFICATION-BASED ASPIRATION LEVEL # (Classification- {">","<","<"})
z_0 = [1.559e-01, 1.561e-01, 9.753e-01];                                    %PARETO OPTIMAL POINT CORRESPONDING TO USER DEFINED START POINT

%% REFERENCE POINTS
dir1 = asp - z_0;
r_p = z_0 + [0.1*dir1;  0.2*dir1; 0.3*dir1; 0.4*dir1; 0.5*dir1];           %IN PR ALGORITHM REFERENCE POINTS ARE PROVIDED ITERATIVELY 

%% CALCULATION OF STEP SIZE 
% # ref_point = z + t*ref_dir
% 
% # SPEED OPTIONS = ["1","2","3","4","5"]
% 
% # SPEED ["5","5","5","5","5"] 
% 
% # t = ["5","5","5","5","5"]*(min(ideal point - nadir point))/(5*10) = [0.1, 0.1, 0.1, 0.1, 0.1]
% 
% # Increement in step size [0.1, 0.1, 0.1, 0.1, 0.1]
% 
% # [t1, t2, t3, t4, t5] = [0.1, 0.2, 0.3, 0.4, 0.5]

%% VALUES OF REFERENCE POINTS
% r_p = [[0.240, 0.141, 0.878]   %1
%        [0.325, 0.125, 0.780]   %2
%        [0.409, 0.109, 0.683]   %3
%        [0.494, 0.094, 0.585]   %4
%        [0.578, 0.078, 0.488]]; %5

%% PREFERRED SOLUTIONS
res_a =[[0.261, 0.153, 0.953] %1 SPEED: `5'                                %IN PR ALGORITHM PREFERRED SOLUTIONS ARE COMPUTED ITERATIVELY
        [0.380, 0.146, 0.913] %2 SPEED: `5'
        [0.509, 0.136, 0.849] %3 SPEED: `5'
        [0.639, 0.123, 0.759] %4 SPEED: `5'
        [0.761, 0.102, 0.641] %5 SPEED: `5'
        ]; 

%% TAKING TURN
asp1 = [0, 1, 0];                                                           %CLASSIFICATION-BASED ASPIRATION LEVEL (Classification- {"<",">","<"})
r_1 = asp1 - res_a(5,:);                                                    %NEW REFERENCE DIRECTION

%% REFERENCE POINTS IN NEW REFERENCE DIRECTION
r_p1 = res_a(5,:) + [r_1*0.1; r_1*0.2; r_1*0.30; r_1*0.36; r_1*0.40];       %REFERENCE POINTS IN NEW REFERENCE DIRECTION
%IN PR ALGORITHM REFERENCE POINTS ARE PROVIDED ITERATIVELY 

%% CALCULATION OF STEP SIZE 
% # ref_point = z + t*ref_dir
% 
% # SPEED OPTIONS = ["1","2","3","4","5"]
% 
% # SPEED ["5","5","5","3","2"] 
% 
% # t = ["5","5","5","3","2"]*(min(ideal point - nadir point))/(5*10) = [0.1, 0.1, 0.1, 0.06, 0.04]
% 
% # Increement in step size [0.1, 0.1, 0.1, 0.06, 0.04]
% 
% # [t1, t2, t3, t4, t5] = [0.1, 0.2, 0.3, 0.36, 0.40]

%% VALUES OF REFERENCE POINTS AFTER TAKING TURN
% r_p1 = [[0.685, 0.192, 0.577]
%         [0.608, 0.282, 0.513]
%         [0.532, 0.371, 0.449]
%         [0.487, 0.425, 0.410]
%         [0.456, 0.461, 0.385]];

%% PREFERRED SOLUTIONS AFETR TAKING TURN
res_b =[[0.747, 0.210, 0.631] %1 SPEED: `5'                                 %IN PR ALGORITHM PREFERRED SOLUTIONS ARE COMPUTED ITERATIVELY 
        [0.721, 0.333, 0.608] %2 SPEED: `5'
        [0.674, 0.471, 0.569] %3 SPEED: `5'
        [0.636, 0.556, 0.536] %4 SPEED: `3'
        [0.605, 0.612, 0.510] %5 SPEED: `2'
        ]; 

%% SCATTER PLOT FOR PREFEREED SOLUTIONS 
scatter3(F(:,1),F(:,2),F(:,3),5,'c','filled','o','LineWidth',2); hold on
scatter3(r_0(:,1),r_0(:,2),r_0(:,3),50,'k','filled','o'); 
scatter3(z_0(:,1),z_0(:,2),z_0(:,3),50,'r','filled','o','MarkerEdgeColor','k','LineWidth',1.5); 
scatter3(r_p(:,1),r_p(:,2),r_p(:,3),30,'b','filled','o','MarkerEdgeColor','b','LineWidth',1);
scatter3(res_a(:,1),res_a(:,2),res_a(:,3),40,'r','filled','o','MarkerEdgeColor','r','LineWidth',1);
scatter3(res_a(5,1),res_a(5,2),res_a(5,3),50,'m','filled','o','MarkerEdgeColor','k','LineWidth',1.5); 
scatter3(res_b(5,1),res_b(5,2),res_b(5,3),40,'g','filled','o','MarkerEdgeColor','k','LineWidth',1.5);
scatter3(r_p1(:,1),r_p1(:,2),r_p1(:,3),30,'b','filled','o','MarkerEdgeColor','b','LineWidth',1);
scatter3(res_b([1 2 3 4],1),res_b([1 2 3 4],2),res_b([1 2 3 4],3),40,'r','filled','o','MarkerEdgeColor','r','LineWidth',1);


legend('Pareto Front', 'Start Point', '$z_0$', 'Ref Points', 'PR solns', 'Turn','Final Solution', 'interpreter','latex','FontSize',12)
xlabel('$f_1$', 'interpreter','latex','FontSize',15); 
ylabel('$f_2$', 'interpreter','latex','FontSize',15);
zlabel('$f_3$', 'interpreter','latex','FontSize',15);
sgtitle('Refer to Figure 1 in manuscript')
