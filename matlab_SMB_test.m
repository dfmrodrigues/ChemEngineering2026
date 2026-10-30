clc
clear
close all
profile on
% u = [8, 3.73447, 4.2561, 3.12978, 0.393381]';
u = [
  800.0000
  372.4387
  424.4895
  312.1199
   39.4438]/100;
mu = [
  509.8402
  196.1648]/100;
mulower = [
    1.9392
         0
         0
         0
         0]/100;
nv = 5;
u = u*ones(1,nv);
% u = [7.24586, 7.39459, 8, 8, 8;...
%     7.24586, 7.39459, 1.48531, 4.1769, 5.25529;...
%     7.24586, 7.39459, 6.82218, 4.17691, 5.25529;...
%     2.92862, 2.2387, 6.82218, 4.17691, 5.25529;...
%     0.276374, 0.276374, 0.276374, 0.276374, 0.276374];
ip = 50;
configuration = [1,1,1,1];
dv = 0;
C_feed = [0.5; 0.5];
K = [0.518; 0.743];
q1 = 0;
q2 = 0;
b1 = [0; 0];
b2 = [0; 0];
k_h = 3.6*6.84/(1-0.389)*[1/0.518; 1/0.743];
L_c = 8/4;
eps = 0.389;
Peinv = 0;
p0 = [C_feed;K;q1;q2;b1;b2;k_h;L_c;eps;Peinv];
idxJ = [1,2,4,5];
rtoly = 1e-3;
atoly = 1e-6;
rtola = 1e-3;
atola = 1e-6;
check = true;
checkp = true;
if check == true || checkp == true
[obj,ineq,gobj,gineq,...
    dobj_du,dineq_du,gpobj,gpineq,dobj_dp,dineq_dp,...
    J_1,J_2,J_3,J_4,J_5,J_6,gJ_1,gJ_2,gJ_3,gJ_4,gJ_5,gJ_6,...
    dJ_1_du,dJ_2_du,dJ_3_du,dJ_4_du,dJ_5_du,dJ_6_du,gpJ_1,gpJ_2,gpJ_3,gpJ_4,gpJ_5,gpJ_6,dJ_1_dp,dJ_2_dp,dJ_3_dp,dJ_4_dp,dJ_5_dp,dJ_6_dp] =...
    matlab_SMB_obj_constr(u,ip,configuration,dv,p0,idxJ,rtoly,atoly,rtola,atola,check,checkp);
else
[obj,ineq,gobj,gineq,...
    ~,~,~,~,~,~,...
    J_1,J_2,J_3,J_4,J_5,J_6,gJ_1,gJ_2,gJ_3,gJ_4,gJ_5,gJ_6,...
    ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] =...
    matlab_SMB_obj_constr(u,ip,configuration,dv,p0,idxJ,rtoly,atoly,rtola,atola);
end
disp(obj)
disp(ineq(1))
disp(ineq(2))
disp(ineq(3))
disp(ineq(4))
disp(ineq(5))
if check == true
disp(gJ_1)
disp(dJ_1_du-gJ_1)
end
if checkp == true
disp(gpJ_1)
disp(dJ_1_dp-gpJ_1)
end
if check == true
disp(gJ_2)
disp(dJ_2_du-gJ_2)
end
if checkp == true
disp(gpJ_2)
disp(dJ_2_dp-gpJ_2)
end
if check == true
disp(gJ_3)
disp(dJ_3_du-gJ_3)
end
if checkp == true
disp(gpJ_3)
disp(dJ_3_dp-gpJ_3)
end
if check == true
disp(gJ_4)
disp(dJ_4_du-gJ_4)
end
if checkp == true
disp(gpJ_4)
disp(dJ_4_dp-gpJ_4)
end
if check == true
disp(gJ_5)
disp(dJ_5_du-gJ_5)
end
if checkp == true
disp(gpJ_5)
disp(dJ_5_dp-gpJ_5)
end
if check == true
disp(gJ_6)
disp(dJ_6_du-gJ_6)
end
if checkp == true
disp(gpJ_6)
disp(dJ_6_dp-gpJ_6)
end
if check == true
disp(gobj)
disp(dobj_du'-gobj)
end
if checkp == true
disp(gpobj)
disp(dobj_dp'-gpobj)
end
if check == true
disp(gineq(:,1))
disp(dineq_du(1,:)'-gineq(:,1))
end
if checkp == true
disp(gpineq(:,1))
disp(dineq_dp(1,:)'-gpineq(:,1))
end
if check == true
disp(gineq(:,2))
disp(dineq_du(2,:)'-gineq(:,2))
end
if checkp == true
disp(gpineq(:,2))
disp(dineq_dp(2,:)'-gpineq(:,2))
end
if check == true
disp(gineq(:,3))
disp(dineq_du(3,:)'-gineq(:,3))
end
if checkp == true
disp(gpineq(:,3))
disp(dineq_dp(3,:)'-gpineq(:,3))
end
if check == true
disp(gineq(:,4))
disp(dineq_du(4,:)'-gineq(:,4))
end
if checkp == true
disp(gpineq(:,4))
disp(dineq_dp(4,:)'-gpineq(:,4))
end
if check == true
disp(gineq(:,5))
disp(dineq_du(5,:)'-gineq(:,5))
end
if checkp == true
disp(gpineq(:,5))
disp(dineq_dp(5,:)'-gpineq(:,5))
end
profile viewer
