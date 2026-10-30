function [obj,ineq,gobj,gineq,...
    dobj_du,dineq_du,gpobj,gpineq,dobj_dp,dineq_dp,...
    J_1,J_2,J_3,J_4,J_5,J_6,gJ_1,gJ_2,gJ_3,gJ_4,gJ_5,gJ_6,...
    dJ_1_du,dJ_2_du,dJ_3_du,dJ_4_du,dJ_5_du,dJ_6_du,gpJ_1,gpJ_2,gpJ_3,gpJ_4,gpJ_5,gpJ_6,dJ_1_dp,dJ_2_dp,dJ_3_dp,dJ_4_dp,dJ_5_dp,dJ_6_dp] =...
    matlab_SMB_obj_constr(u,ip,configuration,dv,p0,idxJ,rtoly,atoly,rtola,atola,check,checkp)
if nargin == 10
    check = false;
    checkp = false;
    if nargout > 4
        fast = true;
    else
        fast = false;
    end
else
    fast = false;
end
nv = size(u,2);
nn = length(ip);
dx = 1./ip;
if nn == 1
N_c = sum(configuration);
else
N_c = sum(configuration)+[0;(dv > 1)+6*(dv > 0)];
end
N_cx = configuration(1);
N_cf = sum(configuration(1:2));
N_cr = sum(configuration(1:3));
n_p1 = floor(1./dx) + 1;
n_p = N_c.*n_p1;
C_feed = p0(1:2);
if ~fast
Cycles = 15*N_c(1);
else
Cycles = 30*N_c(1);
end
t0 = 0;
tf = 1/nv;%N_c(1);
du = 1e-8;
dp = 1e-8;
nJ = 6;
nx = sum(round(1./dx+1).*N_c)*4+nJ;
np = size(p0,1);
nu = size(u,1);
ny = nx+np+nu;
if dv > 1
nyF = n_p1(2)*N_c(1)*4;
end
if dv > 0
nyI = n_p1(2)*5*N_c(1)*4;
end
if ~fast
dt = 0.001;
else
dt = 0.01;
end
reps = ceil(Cycles);
% nt = (tf-t0)/dt;
tspan = t0:dt:tf;
x0 = zeros(nx,1);
y0 = [x0;p0;zeros(nu,1)];
if dv > 1
yF0 = zeros(nyF,1);
end
if dv > 0
yI0 = zeros(nyI,1);
end
if check
y0u = cell(nu,1);
for s = 1:nu
y0u{s} = [x0;p0;zeros(nu,1)]*ones(1,nv);
end
if dv > 1
yF0u = cell(nu,1);
for s = 1:nu
yF0u{s} = zeros(nyF,nv);
end
end
if dv > 0
yI0u = cell(nu,1);
for s = 1:nu
yI0u{s} = zeros(nyI,nv);
end
end
end
if checkp
y0p = [x0*ones(1,np);p0+dp*eye(np);zeros(nu,np)];
if dv > 1
yF0p = zeros(nyF,np);
end
if dv > 0
yI0p = zeros(nyI,np);
end
end
% tout = cell(1,reps);
% yout = cell(1,reps);
tout = cell(nv,reps);
yout = cell(nv,reps);
if check
% youtku = cell(nu,1);
% for s = 1:nu
% youtku{s} = zeros(ny,nt+1,nv);
% end
youtku = cell(nu,nv);
end
if check
% youtkp = zeros(ny,nt+1,np);
youtkp = cell(1,np);
end
% opts = odeset('RelTol',rtoly,'AbsTol',atoly);
for k = 1:reps
    if check
    disp(k)
    end
    % tout{k} = zeros(1,nv*nt+1);
    % yout{k} = zeros(ny,nv*nt+1);
    % tout{k}(1) = t0;
    % yout{k}(:,1) = y0;
    tk = mod((k-1)*nv*(tf-t0),N_c(1));%floor((k-1)/N_c(1))*N_c(1);%(k-1)*(tf-t0);
    if dv > 1
    idx_feed = mod(N_cf-1+tk,N_c(1));
    end
    if dv > 0
    idx_ext = mod(N_cx-1+tk,N_c(1));
    idx_raf = mod(N_cr-1+tk,N_c(1));
    end
    y0(end-nu-np-nJ+1:end-nu-np) = zeros(nJ,1);
    if dv > 1
    y0(n_p(1)*4+n_p1(2)*(N_c(2)-7)*2+(1:(n_p1(2)*2))) = yF0(n_p1(2)*idx_feed*2+(1:(n_p1(2)*2)));
    y0(n_p(1)*4+n_p1(2)*(2*N_c(2)-7)*2+(1:(n_p1(2)*2))) = yF0(n_p1(2)*(N_c(1)+idx_feed)*2+(1:(n_p1(2)*2)));
    end
    if dv > 0
    for j = 1:N_c(1)
        if j-1 == idx_ext || j-1 == idx_raf
            continue
        end
        y0(n_p(1)*4+n_p1(2)*(j-1)*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)));
        y0(n_p(1)*4+n_p1(2)*(N_c(2)+j-1)*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*(N_c(1)+j-1)*2+(1:(n_p1(2)*2)));
    end
    y0(n_p(1)*4+n_p1(2)*(N_c(2)-6)*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*1*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)));
    y0(n_p(1)*4+n_p1(2)*(2*N_c(2)-6)*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*1*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)));
    y0(n_p(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*2*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)));
    y0(n_p(1)*4+n_p1(2)*(N_c(2)+idx_ext)*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*2*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)));
    y0(n_p(1)*4+n_p1(2)*(N_c(2)-3)*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*3*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)));
    y0(n_p(1)*4+n_p1(2)*(2*N_c(2)-3)*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*3*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)));
    y0(n_p(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*4*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)));
    y0(n_p(1)*4+n_p1(2)*(N_c(2)+idx_raf)*2+(1:(n_p1(2)*2))) = yI0(n_p1(2)*4*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)));
    end
    if check
        for s = 1:nu
        y0u{s}(end-nu-np-nJ+1:end-nu-np,:) = zeros(nJ,1)*ones(1,nv);
        if dv > 1
        y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)-7)*2+(1:(n_p1(2)*2)),:) = yF0u{s}(n_p1(2)*idx_feed*2+(1:(n_p1(2)*2)),:);
        y0u{s}(n_p(1)*4+n_p1(2)*(2*N_c(2)-7)*2+(1:(n_p1(2)*2)),:) = yF0u{s}(n_p1(2)*(N_c(1)+idx_feed)*2+(1:(n_p1(2)*2)),:);
        end
        if dv > 0
        for j = 1:N_c(1)
            if j-1 == idx_ext || j-1 == idx_raf
                continue
            end
            y0u{s}(n_p(1)*4+n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:);
            y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)+j-1)*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*(N_c(1)+j-1)*2+(1:(n_p1(2)*2)),:);
        end
        y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)-6)*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*1*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:);
        y0u{s}(n_p(1)*4+n_p1(2)*(2*N_c(2)-6)*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*1*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:);
        y0u{s}(n_p(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*2*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:);
        y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)+idx_ext)*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*2*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:);
        y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)-3)*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*3*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:);
        y0u{s}(n_p(1)*4+n_p1(2)*(2*N_c(2)-3)*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*3*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:);
        y0u{s}(n_p(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*4*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:);
        y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)+idx_raf)*2+(1:(n_p1(2)*2)),:) = yI0u{s}(n_p1(2)*4*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:);
        end
        end
    end
    if checkp
        y0p(end-nu-np-nJ+1:end-nu-np,:) = zeros(nJ,1)*ones(1,np);
        if dv > 1
        y0p(n_p(1)*4+n_p1(2)*(N_c(2)-7)*2+(1:(n_p1(2)*2)),:) = yF0p(n_p1(2)*idx_feed*2+(1:(n_p1(2)*2)),:);
        y0p(n_p(1)*4+n_p1(2)*(2*N_c(2)-7)*2+(1:(n_p1(2)*2)),:) = yF0p(n_p1(2)*(N_c(1)+idx_feed)*2+(1:(n_p1(2)*2)),:);
        end
        if dv > 0
        for j = 1:N_c(1)
            if j-1 == idx_ext || j-1 == idx_raf
                continue
            end
            y0p(n_p(1)*4+n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:);
            y0p(n_p(1)*4+n_p1(2)*(N_c(2)+j-1)*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*(N_c(1)+j-1)*2+(1:(n_p1(2)*2)),:);
        end
        y0p(n_p(1)*4+n_p1(2)*(N_c(2)-6)*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*1*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:);
        y0p(n_p(1)*4+n_p1(2)*(2*N_c(2)-6)*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*1*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:);
        y0p(n_p(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*2*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:);
        y0p(n_p(1)*4+n_p1(2)*(N_c(2)+idx_ext)*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*2*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:);
        y0p(n_p(1)*4+n_p1(2)*(N_c(2)-3)*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*3*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:);
        y0p(n_p(1)*4+n_p1(2)*(2*N_c(2)-3)*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*3*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:);
        y0p(n_p(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*4*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:);
        y0p(n_p(1)*4+n_p1(2)*(N_c(2)+idx_raf)*2+(1:(n_p1(2)*2)),:) = yI0p(n_p1(2)*4*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:);
        end
    end
for i = 1:nv
    tspank = tspan+((k-1)*nv+i-1)*(tf-t0);%tspan+tk;
    y0(end-nu+1:end) = u(:,i);
    if check
        for j = 1:nv
        for s = 1:nu
        y0u{s}(end-nu+1:end,j) = u(:,i);
        y0u{s}(end-nu+s,j) = y0u{s}(end-nu+s,j)+du*(i==j);
        end
        end
    end
    if checkp
        for j = 1:np
        y0p(end-nu+1:end,j) = u(:,i);
        end
    end
    % [toutk,youtk] = ode45f(@matlab_SMB_odefun_y_mex,tspank,y0,opts,tk);
    % toutk = toutk';
    % youtk = youtk';
    % [toutk,youtk] = ode45_rp(@matlab_SMB_odefun_y_mex,tspank,y0,rtoly,atoly,[],tk);
    [toutk,youtk] = ode45_rp(@matlab_SMB_odefun_2_y_mex,[tspank(1),tspank(end)],y0,rtoly,atoly,[],tk);
    % if max(max(abs(yout2k(1:ny,2)-youtk(:,1)))) > 0
    %     keyboard
    % end
    % if max(max(abs(yout2k(end-ny+1:end,end)-youtk(:,end)))) > 0
    %     keyboard
    % end
    if check
    for j = 1:nv
    for s = 1:nu
    % [~,youtku{s}(:,:,j)] = ode45f(@matlab_SMB_odefun_y_mex,tspank,y0u{s}(:,j),opts,tk);
    % youtku{s}(:,:,j) = youtku{s}(:,:,j)';
    % [~,youtku{s}(:,:,j)] = ode45_rp(@matlab_SMB_odefun_y_mex,tspank,y0u{s}(:,j),rtoly,atoly,[],tk);
    [~,youtku{s,j}] = ode45_rp(@matlab_SMB_odefun_2_y_mex,[tspank(1),tspank(end)],y0u{s}(:,j),rtoly,atoly,[],tk);
    end
    end
    end
    if checkp
    for j = 1:np
    % [~,youtkp(:,:,j)] = ode45f(@matlab_SMB_odefun_y_mex,tspank,y0p(:,j),opts,tk);
    % youtkp(:,:,j) = youtkp(:,:,j)';
    % [~,youtkp(:,:,j)] = ode45_rp(@matlab_SMB_odefun_y_mex,tspank,y0p(:,j),rtoly,atoly,[],tk);
    [~,youtkp{j}] = ode45_rp(@matlab_SMB_odefun_2_y_mex,[tspank(1),tspank(end)],y0p(:,j),rtoly,atoly,[],tk);
    end
    end
    % y0 = youtk(:,end);
    y0 = youtk(end-ny+1:end,end);
    if check
    for j = 1:nv
    for s = 1:nu
    % y0u{s}(:,j) = youtku{s}(:,end,j);
    y0u{s}(:,j) = youtku{s,j}(end-ny+1:end,end);
    end
    end
    end
    if checkp
    for j = 1:np
    % y0p(:,j) = youtkp(:,end,j);
    y0p(:,j) = youtkp{j}(end-ny+1:end,end);
    end
    end
    % tout{k}(((i-1)*nt+1):((i)*nt+1)) = toutk;
    % yout{k}(:,((i-1)*nt+1):((i)*nt+1)) = youtk;
    tout{i,k} = toutk;
    yout{i,k} = youtk;
end
    if dv > 1
    yF0(n_p1(2)*idx_feed*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(N_c(2)-7)*2+(1:(n_p1(2)*2)));
    yF0(n_p1(2)*(N_c(1)+idx_feed)*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(2*N_c(2)-7)*2+(1:(n_p1(2)*2)));
    end
    if dv > 0
    for j = 1:N_c(1)
        if j-1 == idx_ext || j-1 == idx_raf
            continue
        end
        yI0(n_p1(2)*(j-1)*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)));
        yI0(n_p1(2)*(N_c(1)+j-1)*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(N_c(2)+j-1)*2+(1:(n_p1(2)*2)));
    end
    yI0(n_p1(2)*1*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(N_c(2)-6)*2+(1:(n_p1(2)*2)));
    yI0(n_p1(2)*1*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(2*N_c(2)-6)*2+(1:(n_p1(2)*2)));
    yI0(n_p1(2)*2*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)));
    yI0(n_p1(2)*2*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(N_c(2)+idx_ext)*2+(1:(n_p1(2)*2)));
    yI0(n_p1(2)*3*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(N_c(2)-3)*2+(1:(n_p1(2)*2)));
    yI0(n_p1(2)*3*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(2*N_c(2)-3)*2+(1:(n_p1(2)*2)));
    yI0(n_p1(2)*4*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)));
    yI0(n_p1(2)*4*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2))) = y0(n_p(1)*4+n_p1(2)*(N_c(2)+idx_raf)*2+(1:(n_p1(2)*2)));
    end
    if check
        for s = 1:nu
        if dv > 1
        yF0u{s}(n_p1(2)*idx_feed*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)-7)*2+(1:(n_p1(2)*2)),:);
        yF0u{s}(n_p1(2)*(N_c(1)+idx_feed)*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(2*N_c(2)-7)*2+(1:(n_p1(2)*2)),:);
        end
        if dv > 0
        for j = 1:N_c(1)
            if j-1 == idx_ext || j-1 == idx_raf
                continue
            end
            yI0u{s}(n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:);
            yI0u{s}(n_p1(2)*(N_c(1)+j-1)*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)+j-1)*2+(1:(n_p1(2)*2)),:);
        end
        yI0u{s}(n_p1(2)*1*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)-6)*2+(1:(n_p1(2)*2)),:);
        yI0u{s}(n_p1(2)*1*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(2*N_c(2)-6)*2+(1:(n_p1(2)*2)),:);
        yI0u{s}(n_p1(2)*2*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:);
        yI0u{s}(n_p1(2)*2*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)+idx_ext)*2+(1:(n_p1(2)*2)),:);
        yI0u{s}(n_p1(2)*3*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)-3)*2+(1:(n_p1(2)*2)),:);
        yI0u{s}(n_p1(2)*3*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(2*N_c(2)-3)*2+(1:(n_p1(2)*2)),:);
        yI0u{s}(n_p1(2)*4*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:);
        yI0u{s}(n_p1(2)*4*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:) = y0u{s}(n_p(1)*4+n_p1(2)*(N_c(2)+idx_raf)*2+(1:(n_p1(2)*2)),:);
        end
        end
    end
    if checkp
        if dv > 1
        yF0p(n_p1(2)*idx_feed*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(N_c(2)-7)*2+(1:(n_p1(2)*2)),:);
        yF0p(n_p1(2)*(N_c(1)+idx_feed)*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(2*N_c(2)-7)*2+(1:(n_p1(2)*2)),:);
        end
        if dv > 0
        for j = 1:N_c(1)
            if j-1 == idx_ext || j-1 == idx_raf
                continue
            end
            yI0p(n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:);
            yI0p(n_p1(2)*(N_c(1)+j-1)*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(N_c(2)+j-1)*2+(1:(n_p1(2)*2)),:);
        end
        yI0p(n_p1(2)*1*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(N_c(2)-6)*2+(1:(n_p1(2)*2)),:);
        yI0p(n_p1(2)*1*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(2*N_c(2)-6)*2+(1:(n_p1(2)*2)),:);
        yI0p(n_p1(2)*2*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:);
        yI0p(n_p1(2)*2*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(N_c(2)+idx_ext)*2+(1:(n_p1(2)*2)),:);
        yI0p(n_p1(2)*3*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(N_c(2)-3)*2+(1:(n_p1(2)*2)),:);
        yI0p(n_p1(2)*3*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(2*N_c(2)-3)*2+(1:(n_p1(2)*2)),:);
        yI0p(n_p1(2)*4*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:);
        yI0p(n_p1(2)*4*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:) = y0p(n_p(1)*4+n_p1(2)*(N_c(2)+idx_raf)*2+(1:(n_p1(2)*2)),:);
        end
    end
end
J_1 = y0(end-nu-np-nJ+1);
if check
J_1u = cell(nu,1);
for s = 1:nu
J_1u{s} = y0u{s}(end-nu-np-nJ+1,:);
end
end
if checkp
J_1p = y0p(end-nu-np-nJ+1,:);
end
J_2 = y0(end-nu-np-nJ+2);
if check
J_2u = cell(nu,1);
for s = 1:nu
J_2u{s} = y0u{s}(end-nu-np-nJ+2,:);
end
end
if checkp
J_2p = y0p(end-nu-np-nJ+2,:);
end
J_3 = y0(end-nu-np-nJ+3);
if check
J_3u = cell(nu,1);
for s = 1:nu
J_3u{s} = y0u{s}(end-nu-np-nJ+3,:);
end
end
if checkp
J_3p = y0p(end-nu-np-nJ+3,:);
end
J_4 = y0(end-nu-np-nJ+4);
if check
J_4u = cell(nu,1);
for s = 1:nu
J_4u{s} = y0u{s}(end-nu-np-nJ+4,:);
end
end
if checkp
J_4p = y0p(end-nu-np-nJ+4,:);
end
J_5 = y0(end-nu-np-nJ+5);
if check
J_5u = cell(nu,1);
for s = 1:nu
J_5u{s} = y0u{s}(end-nu-np-nJ+5,:);
end
end
if checkp
J_5p = y0p(end-nu-np-nJ+5,:);
end
J_6 = y0(end-nu-np-nJ+6);
if check
J_6u = cell(nu,1);
for s = 1:nu
J_6u{s} = y0u{s}(end-nu-np-nJ+6,:);
end
end
if checkp
J_6p = y0p(end-nu-np-nJ+6,:);
end
if check
dJ_1_du = zeros(nu,nv);
for s = 1:nu
dJ_1_du(s,:) = (J_1u{s}-J_1)/du;
end
else
dJ_1_du = [];
end
if checkp
dJ_1_dp = (J_1p-J_1)/dp;
else
dJ_1_dp = [];
end
if check
dJ_2_du = zeros(nu,nv);
for s = 1:nu
dJ_2_du(s,:) = (J_2u{s}-J_2)/du;
end
else
dJ_2_du = [];
end
if checkp
dJ_2_dp = (J_2p-J_2)/dp;
else
dJ_2_dp = [];
end
if check
dJ_3_du = zeros(nu,nv);
for s = 1:nu
dJ_3_du(s,:) = (J_3u{s}-J_3)/du;
end
else
dJ_3_du = [];
end
if checkp
dJ_3_dp = (J_3p-J_3)/dp;
else
dJ_3_dp = [];
end
if check
dJ_4_du = zeros(nu,nv);
for s = 1:nu
dJ_4_du(s,:) = (J_4u{s}-J_4)/du;
end
else
dJ_4_du = [];
end
if checkp
dJ_4_dp = (J_4p-J_4)/dp;
else
dJ_4_dp = [];
end
if check
dJ_5_du = zeros(nu,nv);
for s = 1:nu
dJ_5_du(s,:) = (J_5u{s}-J_5)/du;
end
else
dJ_5_du = [];
end
if checkp
dJ_5_dp = (J_5p-J_5)/dp;
else
dJ_5_dp = [];
end
if check
dJ_6_du = zeros(nu,nv);
for s = 1:nu
dJ_6_du(s,:) = (J_6u{s}-J_6)/du;
end
else
dJ_6_du = [];
end
if checkp
dJ_6_dp = (J_6p-J_6)/dp;
else
dJ_6_dp = [];
end
feed = J_3;
if check
feedu = cell(nu,1);
for s = 1:nu
feedu{s} = J_3u{s};
end
end
if checkp
feedp = J_3p;
end
PUR_x = J_2./(J_1 + J_2);
if check
PUR_xu = cell(nu,1);
for s = 1:nu
PUR_xu{s} = J_2u{s}./(J_1u{s} + J_2u{s});
end
end
if checkp
PUR_xp = J_2p./(J_1p + J_2p);
end
REC_x = J_2./(J_3*C_feed(2));
if check
REC_xu = cell(nu,1);
for s = 1:nu
REC_xu{s} = J_2u{s}./(J_3u{s}*C_feed(2));
end
end
if checkp
REC_xp = J_2p./(J_3p*C_feed(2));
end
eluent = J_6;
if check
eluentu = cell(nu,1);
for s = 1:nu
eluentu{s} = J_6u{s};
end
end
if checkp
eluentp = J_6p;
end
PUR_r = J_4./(J_4 + J_5);
if check
PUR_ru = cell(nu,1);
for s = 1:nu
PUR_ru{s} = J_4u{s}./(J_4u{s} + J_5u{s});
end
end
if checkp
PUR_rp = J_4p./(J_4p + J_5p);
end
REC_r = J_4./(J_3*C_feed(1));
if check
REC_ru = cell(nu,1);
for s = 1:nu
REC_ru{s} = J_4u{s}./(J_3u{s}*C_feed(1));
end
end
if checkp
REC_rp = J_4p./(J_3p*C_feed(1));
end
obj = -feed;
ineq = [-PUR_x;-REC_x;eluent;-PUR_r;-REC_r];
if check
dfeed_du = zeros(nu,nv);
for s = 1:nu
dfeed_du(s,:) = (feedu{s}-feed)/du;
end
end
if checkp
dfeed_dp = (feedp-feed)/dp;
end
if check
dPUR_x_du = zeros(nu,nv);
for s = 1:nu
dPUR_x_du(s,:) = (PUR_xu{s}-PUR_x)/du;
end
end
if checkp
dPUR_x_dp = (PUR_xp-PUR_x)/dp;
end
if check
dREC_x_du = zeros(nu,nv);
for s = 1:nu
dREC_x_du(s,:) = (REC_xu{s}-REC_x)/du;
end
end
if checkp
dREC_x_dp = (REC_xp-REC_x)/dp;
end
if check
deluent_du = zeros(nu,nv);
for s = 1:nu
deluent_du(s,:) = (eluentu{s}-eluent)/du;
end
end
if checkp
deluent_dp = (eluentp-eluent)/dp;
end
if check
dPUR_r_du = zeros(nu,nv);
for s = 1:nu
dPUR_r_du(s,:) = (PUR_ru{s}-PUR_r)/du;
end
end
if checkp
dPUR_r_dp = (PUR_rp-PUR_r)/dp;
end
if check
dREC_r_du = zeros(nu,nv);
for s = 1:nu
dREC_r_du(s,:) = (REC_ru{s}-REC_r)/du;
end
end
if checkp
dREC_r_dp = (REC_rp-REC_r)/dp;
end
if check
dobj_du = -dfeed_du(:)';
dineq_du = [-dPUR_x_du(:)';-dREC_x_du(:)';deluent_du(:)';-dPUR_r_du(:)';-dREC_r_du(:)'];
else
dobj_du = [];
dineq_du = [];
end
if checkp
dobj_dp = -dfeed_dp(:)';
dineq_dp = [-dPUR_x_dp(:)';-dREC_x_dp(:)';deluent_dp(:)';-dPUR_r_dp(:)';-dREC_r_dp(:)'];
else
dobj_dp = [];
dineq_dp = [];
end
a0 = zeros(ny,nJ);
a0(end-nu-np-nJ+1:end-nu-np,:) = eye(nJ);
if dv > 1
aF0 = zeros(nyF,nJ);
end
if dv > 0
aI0 = zeros(nyI,nJ);
end
au = zeros(nu,nJ,nv);
ap = zeros(1,nJ,np);
if ~fast
% aoutk = zeros(ny,nt+1,nJ);
aoutk = cell(1,nJ);
% opts = odeset('RelTol',rtola,'AbsTol',atola);
for k = reps:-1:1
    if check
    disp(k)
    end
    tk = mod((k-1)*nv*(tf-t0),N_c(1));%floor((k-1)/N_c(1))*N_c(1);%(k-1)*(tf-t0);
    if dv > 1
    idx_feed = mod(N_cf-1+tk,N_c(1));
    end
    if dv > 0
    idx_ext = mod(N_cx-1+tk,N_c(1));
    idx_raf = mod(N_cr-1+tk,N_c(1));
    end
    if dv > 1
    a0(n_p(1)*4+n_p1(2)*(N_c(2)-7)*2+(1:(n_p1(2)*2)),:) = aF0(n_p1(2)*idx_feed*2+(1:(n_p1(2)*2)),:);
    a0(n_p(1)*4+n_p1(2)*(2*N_c(2)-7)*2+(1:(n_p1(2)*2)),:) = aF0(n_p1(2)*(N_c(1)+idx_feed)*2+(1:(n_p1(2)*2)),:);
    end
    if dv > 0
    for j = 1:N_c(1)
        if j-1 == idx_ext || j-1 == idx_raf
            continue
        end
        a0(n_p(1)*4+n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:);
        a0(n_p(1)*4+n_p1(2)*(N_c(2)+j-1)*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*(N_c(1)+j-1)*2+(1:(n_p1(2)*2)),:);
    end
    a0(n_p(1)*4+n_p1(2)*(N_c(2)-6)*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*1*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:);
    a0(n_p(1)*4+n_p1(2)*(2*N_c(2)-6)*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*1*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:);
    a0(n_p(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*2*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:);
    a0(n_p(1)*4+n_p1(2)*(N_c(2)+idx_ext)*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*2*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:);
    a0(n_p(1)*4+n_p1(2)*(N_c(2)-3)*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*3*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:);
    a0(n_p(1)*4+n_p1(2)*(2*N_c(2)-3)*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*3*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:);
    a0(n_p(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*4*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:);
    a0(n_p(1)*4+n_p1(2)*(N_c(2)+idx_raf)*2+(1:(n_p1(2)*2)),:) = aI0(n_p1(2)*4*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:);
    end
for i = nv:-1:1
    % tspank = fliplr(tspan+((k-1)*nv+i-1)*(tf-t0));%tspan+tk;
    a0(end-nu+1:end,:) = au(:,:,i);
    for j = idxJ
    % [~,aoutk(:,:,j)] = ode45f(@matlab_SMB_odefun_a_mex,tspank,a0(:,j),opts,tout,yout,tk);
    % aoutk(:,:,j) = aout1k(:,:,j)';
    % [~,aoutk(:,:,j)] = ode45_rp(@matlab_SMB_odefun_a_mex,tspank,a0(:,j),rtola,atola,[],tout{k},yout{k},tk,checkp);
    % a0(:,j) = aoutk(:,end,j);
    [~,aoutk{j}] = ode45_rp(@matlab_SMB_odefun_2_a_mex,fliplr(tout{i,k}),a0(:,j),rtola,atola,[],fliplr(yout{i,k}(:,2:end)),tk,checkp);
    a0(:,j) = aoutk{j}(:,end);
    end
    au(:,:,i) = a0(end-nu+1:end,:);
end
    for j = 1:np
        ap(:,:,j) = a0(end-nu-np+j,:);
    end
    a0(end-nu-np-nJ+1:end-nu-np,:) = zeros(nJ);
    if dv > 1
    aF0(n_p1(2)*idx_feed*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(N_c(2)-7)*2+(1:(n_p1(2)*2)),:);
    aF0(n_p1(2)*(N_c(1)+idx_feed)*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(2*N_c(2)-7)*2+(1:(n_p1(2)*2)),:);
    end
    if dv > 0
    for j = 1:N_c(1)
        if j-1 == idx_ext || j-1 == idx_raf
            continue
        end
        aI0(n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(j-1)*2+(1:(n_p1(2)*2)),:);
        aI0(n_p1(2)*(N_c(1)+j-1)*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(N_c(2)+j-1)*2+(1:(n_p1(2)*2)),:);
    end
    aI0(n_p1(2)*1*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(N_c(2)-6)*2+(1:(n_p1(2)*2)),:);
    aI0(n_p1(2)*1*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(2*N_c(2)-6)*2+(1:(n_p1(2)*2)),:);
    aI0(n_p1(2)*2*N_c(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*idx_ext*2+(1:(n_p1(2)*2)),:);
    aI0(n_p1(2)*2*N_c(1)*4+n_p1(2)*(N_c(1)+idx_ext)*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(N_c(2)+idx_ext)*2+(1:(n_p1(2)*2)),:);
    aI0(n_p1(2)*3*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(N_c(2)-3)*2+(1:(n_p1(2)*2)),:);
    aI0(n_p1(2)*3*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(2*N_c(2)-3)*2+(1:(n_p1(2)*2)),:);
    aI0(n_p1(2)*4*N_c(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*idx_raf*2+(1:(n_p1(2)*2)),:);
    aI0(n_p1(2)*4*N_c(1)*4+n_p1(2)*(N_c(1)+idx_raf)*2+(1:(n_p1(2)*2)),:) = a0(n_p(1)*4+n_p1(2)*(N_c(2)+idx_raf)*2+(1:(n_p1(2)*2)),:);
    end
end
end
au = permute(au,[1,3,2]);
ap = permute(ap,[1,3,2]);
gJ_1 = au(:,:,1);
gJ_2 = au(:,:,2);
gJ_3 = [0;-1/nv;1/nv;0;0]*ones(1,nv);
gJ_4 = au(:,:,4);
gJ_5 = au(:,:,5);
gJ_6 = [1/nv;0;0;-1/nv;0]*ones(1,nv);
gpJ_1 = ap(:,:,1);
gpJ_2 = ap(:,:,2);
gpJ_3 = zeros(1,np);
gpJ_4 = ap(:,:,4);
gpJ_5 = ap(:,:,5);
gpJ_6 = zeros(1,np);
gfeed = gJ_3;
gPUR_x = gJ_2/(J_1 + J_2)-gJ_1*J_2/(J_1 + J_2)^2-gJ_2*J_2/(J_1 + J_2)^2;
gREC_x = gJ_2/(J_3*C_feed(2))-gJ_3*C_feed(2)*J_2/(J_3*C_feed(2))^2;
geluent = gJ_6;
gPUR_r = gJ_4/(J_4 + J_5)-gJ_4*J_4/(J_4 + J_5)^2-gJ_5*J_4/(J_4 + J_5)^2;
gREC_r = gJ_4/(J_3*C_feed(1))-gJ_3*C_feed(1)*J_4/(J_3*C_feed(1))^2;
gpfeed = gpJ_3;
gpPUR_x = gpJ_2/(J_1 + J_2)-gpJ_1*J_2/(J_1 + J_2)^2-gpJ_2*J_2/(J_1 + J_2)^2;
gpREC_x = gpJ_2/(J_3*C_feed(2))-gpJ_3*C_feed(2)*J_2/(J_3*C_feed(2))^2;
gpeluent = gpJ_6;
gpPUR_r = gpJ_4/(J_4 + J_5)-gpJ_4*J_4/(J_4 + J_5)^2-gpJ_5*J_4/(J_4 + J_5)^2;
gpREC_r = gpJ_4/(J_3*C_feed(1))-gpJ_3*C_feed(1)*J_4/(J_3*C_feed(1))^2;
gobj = -gfeed(:);
gineq = [-gPUR_x(:),-gREC_x(:),geluent(:),-gPUR_r(:),-gREC_r(:)];
gpobj = -gpfeed(:);
gpineq = [-gpPUR_x(:),-gpREC_x(:),gpeluent(:),-gpPUR_r(:),-gpREC_r(:)];
end