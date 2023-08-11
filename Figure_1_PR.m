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
r_0 = [0.15, 0.15, 1.25];                                                  %USER DEFINED START POINT
z_0 = [0.123 0.123 0.985];                                                 %PARETO OPTIMAL POINT CORRESPONDING TO USER DEFINED START POINT

%% ASPIRATION LEVEL
asp = [1.0, z_0(2), 1.0];                                                  %CLASSIFICATION-BASED ASPIRATION LEVEL # (Classification- {">","=","<"})

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
% r_p = [[0.211, 0.123, 0.987]   
%        [0.298, 0.123, 0.988]    
%        [0.386, 0.123, 0.99 ]    
%        [0.474, 0.123, 0.991]    
%        [0.561, 0.123, 0.993]];  

%% PREFERRED SOLUTIONS
res_a = [[0.207 0.121 0.971] %1 SPEED: `5''                                 %IN PR ALGORITHM PREFERRED SOLUTIONS ARE COMPUTED ITERATIVELY
         [0.287 0.119 0.951] %2 SPEED: `5'
         [0.361 0.115 0.926] %3 SPEED: `5'
         [0.429 0.111 0.897] %4 SPEED: `5'
         [0.49  0.107 0.865] %5 SPEED: `5                 
         ];

%% TAKING TURN
asp1 = [0, 1, 0];                                                           %CLASSIFICATION-BASED ASPIRATION LEVEL (Classification- {"<",">","<"})
r_1 = asp1 - res_a(end,:);                                                  %NEW REFERENCE DIRECTION

%% REFERENCE POINTS IN NEW REFERENCE DIRECTION
r_p1 = res_a(end,:) + [r_1*0.1; r_1*0.2; r_1*0.28; r_1*0.34; r_1*0.40; r_1*0.46];       %REFERENCE POINTS IN NEW REFERENCE DIRECTION
%IN PR ALGORITHM REFERENCE POINTS ARE PROVIDED ITERATIVELY 

%% CALCULATION OF STEP SIZE 
% # ref_point = z + t*ref_dir
% 
% # SPEED OPTIONS = ["1","2","3","4","5"]
% 
% # SPEED ["5","5","4","3","3","3"] 
% 
% # t = ["5","5","4","3","3","3"]*(min(ideal point - nadir point))/(5*10) = [0.1, 0.1, 0.08, 0.06, 0.06, 0.06]
% 
% # Increement in step size [0.1, 0.1, 0.08, 0.06, 0.06, 0.06]
% 
% # [t1, t2, t3, t4, t5] = [0.1, 0.2, 0.28, 0.34, 0.40, 0.46]

%% VALUES OF REFERENCE POINTS AFTER TAKING TURN
%  r_p1 = [[0.441, 0.196, 0.779] 
%          [0.392, 0.286, 0.692]
%          [0.353, 0.357, 0.623]
%          [0.323, 0.411, 0.571]
%          [0.294, 0.464, 0.519]
%          [0.264, 0.518, 0.467]];

%% PREFERRED SOLUTIONS AFETR TAKING TURN
res_b =[[0.481 0.214 0.850] %1 SPEED: `5'                                   %IN PR ALGORITHM PREFERRED SOLUTIONS ARE COMPUTED ITERATIVELY         
        [0.463 0.338 0.820] %2 SPEED: `5'
        [0.441 0.446 0.779] %3 SPEED: `4'
        [0.418 0.531 0.738] %4 SPEED: `3'
        [0.389 0.614 0.687] %5 SPEED: `3'
        [0.354 0.694 0.627] %6 SPEED: `3'
        ];

%% SCATTER PLOT FOR PREFEREED SOLUTIONS 
scatter3(F(:,1),F(:,2),F(:,3),5,'c','filled','o','LineWidth',2); hold on
scatter3(r_0(:,1),r_0(:,2),r_0(:,3),50,'k','filled','o'); 
scatter3(z_0(:,1),z_0(:,2),z_0(:,3),50,'r','filled','o','MarkerEdgeColor','k','LineWidth',1.5); 
scatter3(r_p(:,1),r_p(:,2),r_p(:,3),30,'b','filled','o','MarkerEdgeColor','b','LineWidth',1);
scatter3(res_a(:,1),res_a(:,2),res_a(:,3),40,'r','filled','o','MarkerEdgeColor','r','LineWidth',1);
scatter3(res_a(end,1),res_a(end,2),res_a(end,3),50,'m','filled','o','MarkerEdgeColor','k','LineWidth',1.5); 
scatter3(res_b(end,1),res_b(end,2),res_b(end,3),40,'g','filled','o','MarkerEdgeColor','k','LineWidth',1.5);
scatter3(r_p1(:,1),r_p1(:,2),r_p1(:,3),30,'b','filled','o','MarkerEdgeColor','b','LineWidth',1);
scatter3(res_b(1:end-1,1),res_b(1:end-1,2),res_b(1:end-1,3),40,'r','filled','o','MarkerEdgeColor','r','LineWidth',1);


legend('Pareto Front', 'Start Point', '$z_0$', 'Ref Points', 'PR solns', 'Turn','Final Solution', 'interpreter','latex','FontSize',12)
xlabel('$f_1$', 'interpreter','latex','FontSize',15); 
ylabel('$f_2$', 'interpreter','latex','FontSize',15);
zlabel('$f_3$', 'interpreter','latex','FontSize',15);
sgtitle('Refer to Figure 1 in manuscript')