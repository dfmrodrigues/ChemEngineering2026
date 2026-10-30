function [u,tau,exitflag,output,nu] = matlab_SMB_fmincon()
last_u = [];
last_obj = [];
last_ineq = [];
last_gobj = [];
last_gineq = [];
nv = 5;
u0 = [8; 5; 5.01; 3.5]*ones(1,nv);
u0 = [u0(:); 1140/60];
% u0 = [u0(:); 1140/3600];
if nv > 1
A = [-1,0,0,1;...
     -1,1,0,0;...
     0,1,-1,0;...
     0,0,-1,1];
else
A = [-1,0,0,1;...
     -1,1,0,0;...
     0,1,-1,0;...
     0,0,-1,1;...
     -1,0,1,0;...
     0,-1,0,1];
end
A = [kron(eye(nv),A),zeros(nv*size(A,1),1)];
if nv > 1
b = [0; 0; 0; 0]*ones(1,nv);
else
b = [0; 0; 0; 0; 0; 0]*ones(1,nv);
end
b = b(:);
% Aeq = [];
% beq = [];
if nv > 1
lb = [0; 0; 0; 0]*ones(1,nv);
else
lb = [8; 0.1; 0.1; 0.1]*ones(1,nv);
end
lb = [lb(:); 0];
ub = [8; 8; 8; 8]*ones(1,nv);
ub = [ub(:); 60];
% ub = [ub(:); 1];
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
%     'OptimalityTolerance',1e-6);
%     % 'SubproblemAlgorithm','cg');
% [u,tau,exitflag,output,nu] = fmincon(@objective,u0,A,b,Aeq,beq,lb,ub,@nonlconstraints,options);
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
idxJ = 1:2;
rtoly = 1e-3;
atoly = 1e-6;
rtola = 1e-3;
atola = 1e-6;
PUR_xmin = 0.9;
REC_xmin = 0.8;
funcs.objective = @objective;
funcs.gradient = @gradient;
funcs.constraints = @constraints;
funcs.jacobian = @jacobian;
funcs.jacobianstructure = @jacobianstructure;
funcs.iterfunc = @iterfunc;
options.lb = lb;
options.ub = ub;
options.cl = -inf*ones(2+size(A,1),1);
options.cu = [zeros(2,1);b];
options.ipopt.print_level = 5;
options.ipopt.tol = 1e-6;
options.ipopt.hessian_approximation = 'limited-memory';
[u, info] = ipopt(u0,funcs,options);
tau = last_obj;
exitflag = info.status;
output = info.eval;
nu.ineqnonlin = info.lambda(1:2);
nu.ineqlin = info.lambda(3:end);
nu.lower = info.zl;
nu.upper = info.zu;
    function obj_constr(u)
        last_u = u;
        u = [reshape(u(1:end-1),[],nv);u(end)*ones(1,nv)];
        % u = reshape(u,[],nv);
        u(end,:) = u(end,:)/60;
        [obj,ineq,gobj,gineq] = matlab_SMB_obj_constr(u,ip,configuration,dv,p0,idxJ,rtoly,atoly,rtola,atola);
        ineq = ineq(1:2)+[PUR_xmin;REC_xmin];
        gineq = gineq(:,1:2);
        gobj(5*(1:nv),:) = gobj(5*(1:nv),:)/60;
        gineq(5*(1:nv),:) = gineq(5*(1:nv),:)/60;
        idx = kron(ones(1,nv),(1:4))+kron(5*(0:(nv-1)),ones(1,4));
        gobj = [gobj(idx,:);sum(gobj(5*(1:nv),:),1)];
        gineq = [gineq(idx,:);sum(gineq(5*(1:nv),:),1)];
        last_obj = obj;
        last_ineq = ineq;
        last_gobj = gobj;
        last_gineq = gineq;
    end
    function [obj,gobj] = objective(u)
        if(~isequal(u,last_u))
            obj_constr(u);
        end
        obj = last_obj;
        gobj = last_gobj;
    end
    function gobj = gradient(u)
        [~,gobj] = objective(u);
    end
    function [ineq,eq,gineq,geq] = nonlconstraints(u)
        if(~isequal(u,last_u))
            obj_constr(u);
        end
        ineq = last_ineq;
        gineq = last_gineq;
        eq = [];
        geq = [];
    end
    function gineq = constraints(u)
        ineq = nonlconstraints(u);
        gineq = [ineq;A*u];
    end
    function gineq = jacobian(u)
        [~,~,gineq] = nonlconstraints(u);
        gineq = sparse([gineq';A]);
    end
    function res = jacobianstructure()
        res = sparse(ones(2+size(A,1),size(A,2)));
    end
    function b = iterfunc(~, ~, ~)
        pause(0);
        b = true;
    end
end
