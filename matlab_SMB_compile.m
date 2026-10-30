% ip - Spaces between column mesh points (-)
% configuration - Columns configuration in the SMB setup (-)
% nc - Number of components (-)
function matlab_SMB_compile(ip, configuration, nc, dv)
clc
close all

%Parameters
nn = length(ip);
dx = 1./ip; %Dimensionless spacing between column mesh points (-)
if nn == 1
N_c = sum(configuration);
else
N_c = sum(configuration)+[0;(dv > 1)+6*(dv > 0)]; %Total number of columns (-)
end
N_cx = configuration(1); %Total number of columns in section 1 (-)
N_cf = sum(configuration(1:2)); %Total number of columns in section 1 and 2 (-)
N_cr = sum(configuration(1:3)); %Total number of columns in sections 1 to 3 (-)
n_p1 = floor(1./dx) + 1; %Number of spatial discretized points in one column (-)
n_p = N_c.*n_p1; %Total number of spatial discretized points in the system (-)

nJ = 6;
nx = sum(round(1./dx+1).*N_c)*4+nJ;
np = 12+3*nn+(dv > 1)+8*(dv > 0);
nu = 5;
ny = nx+np+nu;

% fid = fopen('matlab_SMB_odefun_y.m', 'w');
% 
% fprintf(fid,'function dydt = matlab_SMB_odefun_y(t,y,tk)\n');
% 
% fprintf(fid, 't = t-tk;\n');
% 
% fprintf(fid,'dydt = matlab_SMB_fun_y(t,y);\n');
% 
% fprintf(fid,'end\n');
% 
% fclose(fid);

% fid = fopen('matlab_SMB_fun_y.m', 'w');
% 
% fprintf(fid,'function dydt = matlab_SMB_fun_y(tk,y)\n');

fid = fopen('matlab_SMB_fun_y.c', 'w');

fprintf(fid,'void matlab_SMB_fun_y(double tk, const double y[%d], double dydt[%d]) {\n', ny, ny);

for n = 1:nn
for i = 1:nc
    for j = 1:n_p(n)
        fprintf(fid,'double C_%d_%d_%d;\n', n, i, j);
        fprintf(fid,'double q_%d_%d_%d;\n', n, i, j);
    end
end
end

fprintf(fid,'double J_1;\n');
fprintf(fid,'double J_2;\n');
fprintf(fid,'double J_3;\n');
fprintf(fid,'double J_4;\n');
fprintf(fid,'double J_5;\n');
fprintf(fid,'double J_6;\n');

for i = 1:nc
    fprintf(fid,'double C_feed_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double K_%d;\n', i);
end

fprintf(fid,'double q1;\n');
fprintf(fid,'double q2;\n');

for i = 1:nc
    fprintf(fid,'double b1_%d;\n', i);
end
for i = 1:nc
    fprintf(fid,'double b2_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double k_h_%d;\n', i);
end

for n = 1:nn
fprintf(fid,'double L_c_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double eps_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double Peinv_%d;\n', n);
end

if dv > 1
fprintf(fid,'double Lfeed;\n');
end
if dv > 0
fprintf(fid,'double Lextin;\n');
fprintf(fid,'double Lext1;\n');
fprintf(fid,'double Lext2;\n');
fprintf(fid,'double Lext3;\n');
fprintf(fid,'double Lrafin;\n');
fprintf(fid,'double Lraf1;\n');
fprintf(fid,'double Lraf2;\n');
fprintf(fid,'double Lraf3;\n');
end

for s = 1:4
    fprintf(fid,'double u_%d;\n', s);
end

fprintf(fid,'double t_s;\n\n');

fprintf(fid,'double eps_ratio;\n');

fprintf(fid,'double u_f;\n');

fprintf(fid,'double u_x;\n');

fprintf(fid,'double u_e;\n');

fprintf(fid,'double u_r;\n');

for n = 1:nn
for s = 1:4
    fprintf(fid,'double gamma_i_%d_%d;\n', n, s);
end
end

if dv > 1
fprintf(fid,'double gamma_feed;\n');
end
if dv > 0
fprintf(fid,'double gamma_extin;\n');
fprintf(fid,'double gamma_ext1;\n');
fprintf(fid,'double gamma_ext2;\n');
fprintf(fid,'double gamma_ext3;\n');
fprintf(fid,'double gamma_rafin;\n');
fprintf(fid,'double gamma_raf1;\n');
fprintf(fid,'double gamma_raf2;\n');
fprintf(fid,'double gamma_raf3;\n');
end

for i = 1:nc
    fprintf(fid,'double alpha_%d;\n', i);
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double C_f_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double q_f_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double dCdx_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double dqdx_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

for k = 1:N_c(1)
    for c = 0:2%(ceil(Cycles/N_c)+1)
        tkc = (c-1)*N_c(1)+k;
        fprintf(fid,'int H_%d;\n', tkc + N_c(1));
    end
end

for k = 1:N_c(1)
    fprintf(fid,'int is_section1_%d;\n', k);
    fprintf(fid,'int is_section2_%d;\n', k);
    fprintf(fid,'int is_section3_%d;\n', k);
    fprintf(fid,'int is_section4_%d;\n', k);
end

for k = 1:N_c(1)
    fprintf(fid,'int is_eluent_%d;\n', k);
    fprintf(fid,'int is_connection_%d;\n', k);
    fprintf(fid,'int is_feed_%d;\n', k);
end

for k = 1:N_c(1)
    fprintf(fid,'int is_extract_%d;\n', k);
    fprintf(fid,'int is_raffinate_%d;\n', k);
end

for n = 1:nn
for k = 1:N_c(n)
    fprintf(fid,'double gamma_%d_%d;\n', n, k);
end
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        fprintf(fid,'double C_in_%d_%d_%d;\n', n, i, k);
    end
end
end

for i = 1:nc
    for k = 1:N_c(nn)
        fprintf(fid,'double Cx_out_%d_%d;\n', i, k);
        fprintf(fid,'double Cr_out_%d_%d;\n', i, k);
    end
end

for n = 1:1
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double n1_%d_%d;\n', i, j);
        end
    end
end
for k = 1:N_c(n)
    for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
        fprintf(fid,'double d1_%d;\n', j);
    end
end
end
for n = 1:1
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double n2_%d_%d;\n', i, j);
        end
    end
end
for k = 1:N_c(n)
    for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
        fprintf(fid,'double d2_%d;\n', j);
    end
end
end

for n = 1:1
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double s_%d_%d;\n', i, j);
        end
    end
end
end

fprintf(fid,'double Cam_x;\n');

fprintf(fid,'double Cbm_x;\n');

fprintf(fid,'double Cam_r;\n');

fprintf(fid,'double Cbm_r;\n');

fprintf(fid,'double u_fm;\n');

fprintf(fid,'double u_em;\n');

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        j = ((k-1)*n_p1(n)+1);
        fprintf(fid,'double fL_C_%d_%d_%d;\n', n, i, j);
        fprintf(fid,'double fL_q_%d_%d_%d;\n', n, i, j);
    end
end
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double f_C_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double f_q_%d_%d_%d;\n', n, i, j);
        end
    end
end
for i = 1:nc
    for k = 1:N_c(n)
        j = (k*n_p1(n));
        fprintf(fid,'double fR_C_%d_%d_%d;\n', n, i, j);
        fprintf(fid,'double fR_q_%d_%d_%d;\n', n, i, j);
    end
end
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double s_C_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double s_q_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n))
            fprintf(fid,'double dot_C_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double dot_q_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

fprintf(fid,'double dot_J_1;\n');
fprintf(fid,'double dot_J_2;\n');
fprintf(fid,'double dot_J_3;\n');
fprintf(fid,'double dot_J_4;\n');
fprintf(fid,'double dot_J_5;\n');
fprintf(fid,'double dot_J_6;\n');

for i = 1:nc
    fprintf(fid,'double dot_C_feed_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double dot_K_%d;\n', i);
end

fprintf(fid,'double dot_q1;\n');
fprintf(fid,'double dot_q2;\n');

for i = 1:nc
    fprintf(fid,'double dot_b1_%d;\n', i);
end
for i = 1:nc
    fprintf(fid,'double dot_b2_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double dot_k_h_%d;\n', i);
end

for n = 1:nn
fprintf(fid,'double dot_L_c_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double dot_eps_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double dot_Peinv_%d;\n', n);
end

if dv > 1
fprintf(fid,'double dot_Lfeed;\n');
end
if dv > 0
fprintf(fid,'double dot_Lextin;\n');
fprintf(fid,'double dot_Lext1;\n');
fprintf(fid,'double dot_Lext2;\n');
fprintf(fid,'double dot_Lext3;\n');
fprintf(fid,'double dot_Lrafin;\n');
fprintf(fid,'double dot_Lraf1;\n');
fprintf(fid,'double dot_Lraf2;\n');
fprintf(fid,'double dot_Lraf3;\n');
end

for s = 1:4
    fprintf(fid,'double dot_u_%d;\n', s);
end

fprintf(fid,'double dot_t_s;\n');

%Building the simulation model
    idx = 0;

    for n = 1:nn
    for i = 1:nc
        for j = 1:n_p(n)
            fprintf(fid,'C_%d_%d_%d = y[%d];\n', n, i, j, idx); idx = idx + 1;
            fprintf(fid,'q_%d_%d_%d = y[%d];\n', n, i, j, idx); idx = idx + 1;
        end
    end
    end

    fprintf(fid,'J_1 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_2 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_3 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_4 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_5 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_6 = y[%d];\n', idx); idx = idx + 1;

    for i = 1:nc
        fprintf(fid,'C_feed_%d = y[%d];\n', i, idx); idx = idx + 1;
    end

    for i = 1:nc
        fprintf(fid,'K_%d = y[%d];\n', i, idx); idx = idx + 1;
    end

    fprintf(fid,'q1 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'q2 = y[%d];\n', idx); idx = idx + 1;

    for i = 1:nc
        fprintf(fid,'b1_%d = y[%d];\n', i, idx); idx = idx + 1;
    end
    for i = 1:nc
        fprintf(fid,'b2_%d = y[%d];\n', i, idx); idx = idx + 1;
    end

    for i = 1:nc
        fprintf(fid,'k_h_%d = y[%d];\n', i, idx); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'L_c_%d = y[%d];\n', n, idx); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'eps_%d = y[%d];\n', n, idx); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'Peinv_%d = y[%d];\n', n, idx); idx = idx + 1;
    end

    if dv > 1
    fprintf(fid,'Lfeed = y[%d];\n', idx); idx = idx + 1;
    end
    if dv > 0
    fprintf(fid,'Lextin = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lext1 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lext2 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lext3 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lrafin = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lraf1 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lraf2 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lraf3 = y[%d];\n', idx); idx = idx + 1;
    end

    for s = 1:4
        fprintf(fid,'u_%d = y[%d];\n', s, idx); idx = idx + 1;
    end

    fprintf(fid,'t_s = y[%d];\n\n', idx);

    for n = 1:1
    fprintf(fid,'eps_ratio = (1-eps_%d)/eps_%d;\n', n, n);
    end

    fprintf(fid,'u_f = u_3 - u_2;\n');

    fprintf(fid,'u_x = u_1 - u_2;\n');

    fprintf(fid,'u_e = u_1 - u_4;\n');

    fprintf(fid,'u_r = u_3 - u_4;\n');

    for n = 1:nn
    for s = 1:4
        fprintf(fid,'gamma_i_%d_%d = u_%d*t_s/L_c_%d/eps_%d;\n', n, s, s, n, n);
    end
    end

    if dv > 1
    fprintf(fid,'gamma_feed = u_f*t_s/Lfeed/eps_%d;\n', nn);
    end
    if dv > 0
    fprintf(fid,'gamma_extin = u_1*t_s/Lextin/eps_%d;\n', nn);
    fprintf(fid,'gamma_ext1 = u_1*t_s/Lext1/eps_%d;\n', nn);
    fprintf(fid,'gamma_ext2 = u_2*t_s/Lext2/eps_%d;\n', nn);
    fprintf(fid,'gamma_ext3 = u_2*t_s/Lext3/eps_%d;\n', nn);
    fprintf(fid,'gamma_rafin = u_3*t_s/Lrafin/eps_%d;\n', nn);
    fprintf(fid,'gamma_raf1 = u_3*t_s/Lraf1/eps_%d;\n', nn);
    fprintf(fid,'gamma_raf2 = u_4*t_s/Lraf2/eps_%d;\n', nn);
    fprintf(fid,'gamma_raf3 = u_4*t_s/Lraf3/eps_%d;\n', nn);
    end

    for i = 1:nc
        fprintf(fid,'alpha_%d = k_h_%d*t_s;\n', i, i);
    end

    %Declaring the derivatives and interpolations
    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                jf = j + 1;
                fprintf(fid,'C_f_%d_%d_%d = (C_%d_%d_%d + C_%d_%d_%d)*0.5;\n', n, i, j, n, i, jf, n, i, j);
                fprintf(fid,'q_f_%d_%d_%d = (q_%d_%d_%d + q_%d_%d_%d)*0.5;\n', n, i, j, n, i, jf, n, i, j);
            end
        end
    end
    end

    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                jf = j + 1;
                fprintf(fid,'dCdx_%d_%d_%d = (C_%d_%d_%d - C_%d_%d_%d)*%d;\n', n, i, j, n, i, jf, n, i, j, ip(n));
                fprintf(fid,'dqdx_%d_%d_%d = (q_%d_%d_%d - q_%d_%d_%d)*%d;\n', n, i, j, n, i, jf, n, i, j, ip(n));
            end
        end
    end
    end

    for k = 1:N_c(1)
        for c = 0:2%(ceil(Cycles/N_c)+1)
            tkc = (c-1)*N_c(1)+k;
            if tkc >= 0
                fprintf(fid,'H_%d = (tk-%d) >= 0;\n', tkc + N_c(1), tkc);
    %             fprintf(fid,'H_%d = 1/(1+exp(-w*(t-%d)));\n', tkc + N_c(1), tkc);
            else
                fprintf(fid,'H_%d = (tk+%d) >= 0;\n', tkc + N_c(1), -tkc);
    %             fprintf(fid,'H_%d = 1/(1+exp(-w*(t+%d)));\n', tkc + N_c(1), -tkc);
            end
        end
    end

    for k = 1:N_c(1)
        fprintf(fid,'is_section1_%d = (1 - H_%d + H_%d - H_%d + H_%d);\n', k, k, k + N_c(1) - N_cx, k + N_c(1), k + N_c(1) + N_c(1) - N_cx);
        fprintf(fid,'is_section2_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cf, k + N_c(1) - N_cx, k + N_c(1) + N_c(1) - N_cf, k + N_c(1) + N_c(1) - N_cx);
        fprintf(fid,'is_section3_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cr, k + N_c(1) - N_cf, k + N_c(1) + N_c(1) - N_cr, k + N_c(1) + N_c(1) - N_cf);
        fprintf(fid,'is_section4_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k, k + N_c(1) - N_cr, k + N_c(1), k + N_c(1) + N_c(1) - N_cr);
    end
    
    for k = 1:N_c(1)
        fprintf(fid,'is_eluent_%d = (1 - H_%d + H_%d - H_%d + H_%d);\n', k, k, k + N_c(1) - 1, k + N_c(1), k + N_c(1) + N_c(1) - 1);
        fprintf(fid,'is_connection_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k, k + N_c(1) - N_cf - 1, k + N_c(1), k + N_c(1) + N_c(1) - N_cf - 1);
        fprintf(fid,'is_connection_%d = is_connection_%d + (H_%d - H_%d + H_%d - H_%d);\n', k, k, k + N_c(1) - N_cf, k + N_c(1) - 1, k + N_c(1) + N_c(1) - N_cf, k + N_c(1) + N_c(1) - 1);
        fprintf(fid,'is_feed_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cf - 1, k + N_c(1) - N_cf, k + N_c(1) + N_c(1) - N_cf - 1, k + N_c(1) + N_c(1) - N_cf);
    end
    
    %Quadrature functions terms
    for k = 1:N_c(1)
        fprintf(fid,'is_extract_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cx - 1, k + N_c(1) - N_cx, k + N_c(1) + N_c(1) - N_cx - 1, k + N_c(1) + N_c(1) - N_cx);
        fprintf(fid,'is_raffinate_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cr - 1, k + N_c(1) - N_cr, k + N_c(1) + N_c(1) - N_cr - 1, k + N_c(1) + N_c(1) - N_cr);
    end

    for n = 1:nn
    for k = 1:N_c(1)
        fprintf(fid,'gamma_%d_%d = 0;\n', n, k);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_section1_%d*gamma_i_%d_1;\n', n, k, n, k, k, n);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_section2_%d*gamma_i_%d_2;\n', n, k, n, k, k, n);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_section3_%d*gamma_i_%d_3;\n', n, k, n, k, k, n);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_section4_%d*gamma_i_%d_4;\n', n, k, n, k, k, n);
        if n == nn
        if dv > 0
        if k == N_c(1)
            kf = 1;
        else
            kf = k + 1;
        end
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_extract_%d*(gamma_ext3-gamma_i_%d_1);\n', n, k, n, k, kf, n);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_raffinate_%d*(gamma_raf3-gamma_i_%d_3);\n', n, k, n, k, kf, n);
        end
        end
    end
    if n == nn
    if dv > 1
    fprintf(fid,'gamma_%d_%d = gamma_feed;\n', n, N_c(nn)-6);
    end
    if dv > 0
    fprintf(fid,'gamma_%d_%d = gamma_extin;\n', n, N_c(nn)-5);
    fprintf(fid,'gamma_%d_%d = gamma_ext1;\n', n, N_c(nn)-4);
    fprintf(fid,'gamma_%d_%d = gamma_ext2;\n', n, N_c(nn)-3);
    fprintf(fid,'gamma_%d_%d = gamma_rafin;\n', n, N_c(nn)-2);
    fprintf(fid,'gamma_%d_%d = gamma_raf1;\n', n, N_c(nn)-1);
    fprintf(fid,'gamma_%d_%d = gamma_raf2;\n', n, N_c(nn));
    end
    end
    end
    
    for n = 1:nn
    nf = nn + 1 - n;
    for i = 1:nc
        for k = 1:N_c(nf)
            if n == nn
            if k == 1
                j = N_c(1)*n_p1(n);
            else
                j = (k-1)*n_p1(n);
            end
            fprintf(fid,'C_in_%d_%d_%d = 0;\n', nf, i, k);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_eluent_%d*C_%d_%d_%d*u_4/u_1;\n', nf, i, k, nf, i, k, k, n, i, j);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_connection_%d*C_%d_%d_%d;\n', nf, i, k, nf, i, k, k, n, i, j);
            if dv > 1
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_feed_%d*(u_f*C_%d_%d_%d/u_3 + u_2*C_%d_%d_%d/u_3);\n', nf, i, k, nf, i, k, k, n, i, (N_c(nn)-6)*n_p1(n), n, i, j);
            else
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_feed_%d*(u_f*C_feed_%d/u_3 + u_2*C_%d_%d_%d/u_3);\n', nf, i, k, nf, i, k, k, i, n, i, j);
            end
            else
            if dv > 1
            if k == N_c(nn)-6
            fprintf(fid,'C_in_%d_%d_%d = C_feed_%d;\n', nf, i, k, i);
            end
            end
            if dv > 0
            if k <= N_c(1)
            if k == N_c(1)
                kf = 1;
            else
                kf = k + 1;
            end
            j = k*n_p1(n);
            fprintf(fid,'C_in_%d_%d_%d = (1-is_extract_%d-is_raffinate_%d)*C_%d_%d_%d;\n', nf, i, k, kf, kf, n, i, j);
            elseif k ~= N_c(nn)-6 && k ~= N_c(nn)-5 && k ~= N_c(nn)-2
            j = (k-1)*n_p1(nf);
            fprintf(fid,'C_in_%d_%d_%d = C_%d_%d_%d;\n', nf, i, k, nf, i, j);
            elseif k ~= N_c(nn)-6
            fprintf(fid,'C_in_%d_%d_%d = 0;\n', nf, i, k);
            end
            else
            j = k*n_p1(n);
            fprintf(fid,'C_in_%d_%d_%d = C_%d_%d_%d;\n', nf, i, k, n, i, j);
            end
            end
        end
        for k = 1:N_c(nf)
            if n == nn
            if dv > 0
            if k == N_c(1)
                kf = 1;
            else
                kf = k + 1;
            end
            j = k*n_p1(nf);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_extract_%d*C_%d_%d_%d;\n', n, i, N_c(nn)-5, n, i, N_c(nn)-5, kf, nf, i, j);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_raffinate_%d*C_%d_%d_%d;\n', n, i, N_c(nn)-2, n, i, N_c(nn)-2, kf, nf, i, j);
            j = (N_c(nn)-3)*n_p1(n);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_extract_%d*C_%d_%d_%d;\n', n, i, k, n, i, k, kf, n, i, j);
            j = (N_c(nn))*n_p1(n);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_raffinate_%d*C_%d_%d_%d;\n', n, i, k, n, i, k, kf, n, i, j);
            end
            end
        end
    end
    end
    
    %Quadrature functions terms
    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            if n == nn
            if dv > 0
            if k == N_c(nn)-3
            j = (k-1)*n_p1(n);
            fprintf(fid,'Cx_out_%d_%d = C_%d_%d_%d;\n', i, k, n, i, j);
            else
            fprintf(fid,'Cx_out_%d_%d = 0;\n', i, k);
            end
            if k == N_c(nn)
            j = (k-1)*n_p1(n);
            fprintf(fid,'Cr_out_%d_%d = C_%d_%d_%d;\n', i, k, n, i, j);
            else
            fprintf(fid,'Cr_out_%d_%d = 0;\n', i, k);
            end
            else
            if k <= N_c(1)
            if k == 1
                j = N_c(1)*n_p1(n);
            else
                j = (k-1)*n_p1(n);
            end
            fprintf(fid,'Cx_out_%d_%d = is_extract_%d*C_%d_%d_%d;\n', i, k, k, n, i, j);
            fprintf(fid,'Cr_out_%d_%d = is_raffinate_%d*C_%d_%d_%d;\n', i, k, k, n, i, j);
            end
            end
            end
        end
    end
    end

    fprintf(fid,'if(q1 > 0.0) {\n');
    for n = 1:1
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'n1_%d_%d = b1_%d*C_f_%d_%d_%d;\n', i, j, i, n, i, j);
            end
        end
    end
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'d1_%d = 1+b1_1*C_f_%d_1_%d+b1_2*C_f_%d_2_%d;\n', j, n, j, n, j);
        end
    end
    end
    fprintf(fid,'}\n');
    fprintf(fid,'if(q2 > 0.0) {\n');
    for n = 1:1
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'n2_%d_%d = b2_%d*C_f_%d_%d_%d;\n', i, j, i, n, i, j);
            end
        end
    end
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'d2_%d = 1+b2_1*C_f_%d_1_%d+b2_2*C_f_%d_2_%d;\n', j, n, j, n, j);
        end
    end
    end
    fprintf(fid,'}\n');

    %Defining the isotherms
    for n = 1:1
    for i = 1:nc
        fprintf(fid,'if(K_%d > 0.0) {\n', i);
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'s_%d_%d = alpha_%d*(K_%d*C_f_%d_%d_%d - q_f_%d_%d_%d);\n', i, j, i, i, n, i, j, n, i, j);
            end
        end
        fprintf(fid,'} else {\n');
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'s_%d_%d = -alpha_%d*q_f_%d_%d_%d;\n', i, j, i, n, i, j);
            end
        end
        fprintf(fid,'}\n');
    end
    fprintf(fid,'if(q1 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'s_%d_%d += alpha_%d*q1*n1_%d_%d/d1_%d;\n', i, j, i, i, j, j);
            end
        end
    end
    fprintf(fid,'}\n');
    fprintf(fid,'if(q2 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'s_%d_%d += alpha_%d*q2*n2_%d_%d/d2_%d;\n', i, j, i, i, j, j);
            end
        end
    end
    fprintf(fid,'}\n');
    end

    fprintf(fid,'Cam_x = u_x*(');
    for k = 1:N_c(nn)%(1)
        if k < N_c(nn)%(1)
            fprintf(fid,'Cx_out_1_%d + ',k);
        else
            fprintf(fid,'Cx_out_1_%d',k);
        end
    end
    fprintf(fid,');\n');

    fprintf(fid,'Cbm_x = u_x*(');
    for k = 1:N_c(nn)%(1)
        if k < N_c(nn)%(1)
            fprintf(fid,'Cx_out_2_%d + ',k);
        else
            fprintf(fid,'Cx_out_2_%d',k);
        end
    end
    fprintf(fid,');\n');

    fprintf(fid,'Cam_r = u_r*(');
    for k = 1:N_c(nn)%(1)
        if k < N_c(nn)%(1)
            fprintf(fid,'Cr_out_1_%d + ',k);
        else
            fprintf(fid,'Cr_out_1_%d',k);
        end
    end
    fprintf(fid,');\n');

    fprintf(fid,'Cbm_r = u_r*(');
    for k = 1:N_c(nn)%(1)
        if k < N_c(nn)%(1)
            fprintf(fid,'Cr_out_2_%d + ',k);
        else
            fprintf(fid,'Cr_out_2_%d',k);
        end
    end
    fprintf(fid,');\n');

    fprintf(fid,'u_fm = u_f;\n');

    fprintf(fid,'u_em = u_e;\n');

    %Defining the flux and source terms
    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            j = ((k-1)*n_p1(n)+1);
            fprintf(fid,'fL_C_%d_%d_%d = gamma_%d_%d*C_in_%d_%d_%d;\n', n, i, j, n, k, n, i, k);
            fprintf(fid,'fL_q_%d_%d_%d = 0;\n', n, i, j);
        end
    end
    fprintf(fid,'if(Peinv_%d > 0.0) {\n', n);
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'f_C_%d_%d_%d = gamma_%d_%d*(C_f_%d_%d_%d-dCdx_%d_%d_%d*Peinv_%d);\n', n, i, j, n, k, n, i, j, n, i, j, n);
                fprintf(fid,'f_q_%d_%d_%d = 0;\n', n, i, j);
            end
        end
    end
    fprintf(fid,'} else {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'f_C_%d_%d_%d = gamma_%d_%d*C_f_%d_%d_%d;\n', n, i, j, n, k, n, i, j);
                fprintf(fid,'f_q_%d_%d_%d = 0;\n', n, i, j);
            end
        end
    end
    fprintf(fid,'}\n');
    for i = 1:nc
        for k = 1:N_c(n)
            j = (k*n_p1(n));
            fprintf(fid,'fR_C_%d_%d_%d = gamma_%d_%d*C_%d_%d_%d;\n', n, i, j, n, k, n, i, j);
            fprintf(fid,'fR_q_%d_%d_%d = 0;\n', n, i, j);
        end
    end
    end
    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                if n == 1
                fprintf(fid,'s_C_%d_%d_%d = -eps_ratio*s_%d_%d;\n', n, i, j, i, j);
                fprintf(fid,'s_q_%d_%d_%d = s_%d_%d;\n', n, i, j, i, j);
                else
                fprintf(fid,'s_C_%d_%d_%d = 0;\n', n, i, j);
                fprintf(fid,'s_q_%d_%d_%d = 0;\n', n, i, j);
                end
            end
        end
    end
    end

    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n))
                if j == ((k-1)*n_p1(n)+1)
                    %Danckwerts boundary condition (x = 0)
                    fprintf(fid,'dot_C_%d_%d_%d = (fL_C_%d_%d_%d - f_C_%d_%d_%d)*%d + s_C_%d_%d_%d;\n', n, i, j, n, i, j, n, i, j, 2*ip(n), n, i, j);
                    fprintf(fid,'dot_q_%d_%d_%d = (fL_q_%d_%d_%d - f_q_%d_%d_%d)*%d + s_q_%d_%d_%d;\n', n, i, j, n, i, j, n, i, j, 2*ip(n), n, i, j);
                elseif j == (k*n_p1(n))
                    %Boundary condition at x = 1 (dC/dx = 0)
                    fprintf(fid,'dot_C_%d_%d_%d = (f_C_%d_%d_%d - fR_C_%d_%d_%d)*%d + s_C_%d_%d_%d;\n', n, i, j, n, i, j - 1, n, i, j, 2*ip(n), n, i, j - 1);
                    fprintf(fid,'dot_q_%d_%d_%d = (f_q_%d_%d_%d - fR_q_%d_%d_%d)*%d + s_q_%d_%d_%d;\n', n, i, j, n, i, j - 1, n, i, j, 2*ip(n), n, i, j - 1);
                else
                    %Dynamic equation of C and q at 0 < x < 1
                    fprintf(fid,'dot_C_%d_%d_%d = (f_C_%d_%d_%d - f_C_%d_%d_%d)*%d + 0.5*(s_C_%d_%d_%d + s_C_%d_%d_%d);\n', n, i, j, n, i, j - 1, n, i, j, ip(n), n, i, j, n, i, j - 1);
                    fprintf(fid,'dot_q_%d_%d_%d = (f_q_%d_%d_%d - f_q_%d_%d_%d)*%d + 0.5*(s_q_%d_%d_%d + s_q_%d_%d_%d);\n', n, i, j, n, i, j - 1, n, i, j, ip(n), n, i, j, n, i, j - 1);
                end
            end
        end
    end
    end

    %Quadrature functions
    fprintf(fid,'dot_J_1 = Cam_x;\n');
    fprintf(fid,'dot_J_2 = Cbm_x;\n');
    fprintf(fid,'dot_J_3 = u_fm;\n');
    fprintf(fid,'dot_J_4 = Cam_r;\n');
    fprintf(fid,'dot_J_5 = Cbm_r;\n');
    fprintf(fid,'dot_J_6 = u_em;\n');

    for i = 1:nc
        fprintf(fid,'dot_C_feed_%d = 0;\n', i);
    end

    for i = 1:nc
        fprintf(fid,'dot_K_%d = 0;\n', i);
    end

    fprintf(fid,'dot_q1 = 0;\n');
    fprintf(fid,'dot_q2 = 0;\n');

    for i = 1:nc
        fprintf(fid,'dot_b1_%d = 0;\n', i);
    end
    for i = 1:nc
        fprintf(fid,'dot_b2_%d = 0;\n', i);
    end

    for i = 1:nc
        fprintf(fid,'dot_k_h_%d = 0;\n', i);
    end

    for n = 1:nn
    fprintf(fid,'dot_L_c_%d = 0;\n', n);
    end

    for n = 1:nn
    fprintf(fid,'dot_eps_%d = 0;\n', n);
    end

    for n = 1:nn
    fprintf(fid,'dot_Peinv_%d = 0;\n', n);
    end

    if dv > 1
    fprintf(fid,'dot_Lfeed = 0;\n');
    end
    if dv > 0
    fprintf(fid,'dot_Lextin = 0;\n');
    fprintf(fid,'dot_Lext1 = 0;\n');
    fprintf(fid,'dot_Lext2 = 0;\n');
    fprintf(fid,'dot_Lext3 = 0;\n');
    fprintf(fid,'dot_Lrafin = 0;\n');
    fprintf(fid,'dot_Lraf1 = 0;\n');
    fprintf(fid,'dot_Lraf2 = 0;\n');
    fprintf(fid,'dot_Lraf3 = 0;\n');
    end

    for s = 1:4
        fprintf(fid,'dot_u_%d = 0;\n', s);
    end

    fprintf(fid,'dot_t_s = 0;\n');

    idx = 0;

    % fprintf(fid,'dydt = zeros(size(y));\n');

    % fprintf(fid,'dydt[%d] = dot_t;\n', idx); idx = idx + 1;
    for n = 1:nn
    for i = 1:nc
        for j = 1:n_p(n)
            fprintf(fid,'dydt[%d] = dot_C_%d_%d_%d;\n', idx, n, i, j); idx = idx + 1;
            fprintf(fid,'dydt[%d] = dot_q_%d_%d_%d;\n', idx, n, i, j); idx = idx + 1;
        end
    end
    end

    fprintf(fid,'dydt[%d] = dot_J_1;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_J_2;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_J_3;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_J_4;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_J_5;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_J_6;\n', idx); idx = idx + 1;

    for i = 1:nc
        fprintf(fid,'dydt[%d] = dot_C_feed_%d;\n', idx, i); idx = idx + 1;
    end

    for i = 1:nc
        fprintf(fid,'dydt[%d] = dot_K_%d;\n', idx, i); idx = idx + 1;
    end

    fprintf(fid,'dydt[%d] = dot_q1;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_q2;\n', idx); idx = idx + 1;

    for i = 1:nc
        fprintf(fid,'dydt[%d] = dot_b1_%d;\n', idx, i); idx = idx + 1;
    end
    for i = 1:nc
        fprintf(fid,'dydt[%d] = dot_b2_%d;\n', idx, i); idx = idx + 1;
    end

    for i = 1:nc
        fprintf(fid,'dydt[%d] = dot_k_h_%d;\n', idx, i); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'dydt[%d] = dot_L_c_%d;\n', idx, n); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'dydt[%d] = dot_eps_%d;\n', idx, n); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'dydt[%d] = dot_Peinv_%d;\n', idx, n); idx = idx + 1;
    end

    if dv > 1
    fprintf(fid,'dydt[%d] = dot_Lfeed;\n', idx); idx = idx + 1;
    end
    if dv > 0
    fprintf(fid,'dydt[%d] = dot_Lextin;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_Lext1;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_Lext2;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_Lext3;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_Lrafin;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_Lraf1;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_Lraf2;\n', idx); idx = idx + 1;
    fprintf(fid,'dydt[%d] = dot_Lraf3;\n', idx); idx = idx + 1;
    end

    for s = 1:4
        fprintf(fid,'dydt[%d] = dot_u_%d;\n', idx, s); idx = idx + 1;
    end

    fprintf(fid,'dydt[%d] = dot_t_s;\n', idx);

% fprintf(fid,'end\n');

fprintf(fid,'}\n');

fclose(fid);

% fid = fopen('matlab_SMB_odefun_a.m', 'w');
% 
% fprintf(fid,'function dadt = matlab_SMB_odefun_a(t,a,t_int,y_int,tk)\n');
% 
% fprintf(fid, 'y = binterp(t_int,y_int,t); t = t-tk;\n');
% 
% fprintf(fid,'dadt = matlab_SMB_fun_a(t,a,y);\n');
% 
% fprintf(fid,'end\n');
% 
% fclose(fid);

% fid = fopen('matlab_SMB_fun_a.m', 'w');
% 
% fprintf(fid,'function dadt = matlab_SMB_fun_a(tk,a,y)\n');

fid = fopen('matlab_SMB_fun_a.c', 'w');

fprintf(fid,'void matlab_SMB_fun_a(double tk, const double a[%d], const double y[%d], double dadt[%d], int checkp) {\n', ny, ny, ny);

for n = 1:nn
for i = 1:nc
    for j = 1:n_p(n)
        fprintf(fid,'double C_%d_%d_%d;\n', n, i, j);
        fprintf(fid,'double q_%d_%d_%d;\n', n, i, j);
    end
end
end

fprintf(fid,'double J_1;\n');
fprintf(fid,'double J_2;\n');
fprintf(fid,'double J_3;\n');
fprintf(fid,'double J_4;\n');
fprintf(fid,'double J_5;\n');
fprintf(fid,'double J_6;\n');

for i = 1:nc
    fprintf(fid,'double C_feed_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double K_%d;\n', i);
end

fprintf(fid,'double q1;\n');
fprintf(fid,'double q2;\n');

for i = 1:nc
    fprintf(fid,'double b1_%d;\n', i);
end
for i = 1:nc
    fprintf(fid,'double b2_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double k_h_%d;\n', i);
end

for n = 1:nn
fprintf(fid,'double L_c_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double eps_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double Peinv_%d;\n', n);
end

if dv > 1
fprintf(fid,'double Lfeed;\n');
end
if dv > 0
fprintf(fid,'double Lextin;\n');
fprintf(fid,'double Lext1;\n');
fprintf(fid,'double Lext2;\n');
fprintf(fid,'double Lext3;\n');
fprintf(fid,'double Lrafin;\n');
fprintf(fid,'double Lraf1;\n');
fprintf(fid,'double Lraf2;\n');
fprintf(fid,'double Lraf3;\n');
end

for s = 1:4
    fprintf(fid,'double u_%d;\n', s);
end

fprintf(fid,'double t_s;\n\n');

fprintf(fid,'double eps_ratio;\n');

fprintf(fid,'double u_f;\n');

fprintf(fid,'double u_x;\n');

fprintf(fid,'double u_e;\n');

fprintf(fid,'double u_r;\n');

for n = 1:nn
for s = 1:4
    fprintf(fid,'double gamma_i_%d_%d;\n', n, s);
end
end

if dv > 1
fprintf(fid,'double gamma_feed;\n');
end
if dv > 0
fprintf(fid,'double gamma_extin;\n');
fprintf(fid,'double gamma_ext1;\n');
fprintf(fid,'double gamma_ext2;\n');
fprintf(fid,'double gamma_ext3;\n');
fprintf(fid,'double gamma_rafin;\n');
fprintf(fid,'double gamma_raf1;\n');
fprintf(fid,'double gamma_raf2;\n');
fprintf(fid,'double gamma_raf3;\n');
end

for i = 1:nc
    fprintf(fid,'double alpha_%d;\n', i);
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double C_f_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double q_f_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double dCdx_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double dqdx_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

for k = 1:N_c(1)
    for c = 0:2%(ceil(Cycles/N_c)+1)
        tkc = (c-1)*N_c(1)+k;
        fprintf(fid,'int H_%d;\n', tkc + N_c(1));
    end
end

for k = 1:N_c(1)
    fprintf(fid,'int is_section1_%d;\n', k);
    fprintf(fid,'int is_section2_%d;\n', k);
    fprintf(fid,'int is_section3_%d;\n', k);
    fprintf(fid,'int is_section4_%d;\n', k);
end

for k = 1:N_c(1)
    fprintf(fid,'int is_eluent_%d;\n', k);
    fprintf(fid,'int is_connection_%d;\n', k);
    fprintf(fid,'int is_feed_%d;\n', k);
end

for k = 1:N_c(1)
    fprintf(fid,'int is_extract_%d;\n', k);
    fprintf(fid,'int is_raffinate_%d;\n', k);
end

for n = 1:nn
for k = 1:N_c(n)
    fprintf(fid,'double gamma_%d_%d;\n', n, k);
end
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        fprintf(fid,'double C_in_%d_%d_%d;\n', n, i, k);
    end
end
end

for i = 1:nc
    for k = 1:N_c(nn)
        fprintf(fid,'double Cx_out_%d_%d;\n', i, k);
        fprintf(fid,'double Cr_out_%d_%d;\n', i, k);
    end
end

for n = 1:1
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double n1_%d_%d;\n', i, j);
        end
    end
end
for k = 1:N_c(n)
    for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
        fprintf(fid,'double d1_%d;\n', j);
    end
end
end
for n = 1:1
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double n2_%d_%d;\n', i, j);
        end
    end
end
for k = 1:N_c(n)
    for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
        fprintf(fid,'double d2_%d;\n', j);
    end
end
end

for n = 1:1
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double s_%d_%d;\n', i, j);
        end
    end
end
end

for n = 1:nn
for i = 1:nc
    for j = 1:n_p(n)
        fprintf(fid,'double dT_dC_%d_%d_%d;\n', n, i, j);
        fprintf(fid,'double dT_dq_%d_%d_%d;\n', n, i, j);
    end
end
end

fprintf(fid,'double dT_dJ_1;\n');
fprintf(fid,'double dT_dJ_2;\n');
fprintf(fid,'double dT_dJ_3;\n');
fprintf(fid,'double dT_dJ_4;\n');
fprintf(fid,'double dT_dJ_5;\n');
fprintf(fid,'double dT_dJ_6;\n');

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        j = ((k-1)*n_p1(n)+1);
        fprintf(fid,'double ddot_T_dfL_C_%d_%d_%d;\n', n, i, j);
        fprintf(fid,'double ddot_T_dfL_q_%d_%d_%d;\n', n, i, j);
    end
end
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double ddot_T_df_C_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double ddot_T_df_q_%d_%d_%d;\n', n, i, j);
        end
    end
end
for i = 1:nc
    for k = 1:N_c(n)
        j = (k*n_p1(n));
        fprintf(fid,'double ddot_T_dfR_C_%d_%d_%d;\n', n, i, j);
        fprintf(fid,'double ddot_T_dfR_q_%d_%d_%d;\n', n, i, j);
    end
end
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double ddot_T_ds_C_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double ddot_T_ds_q_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

fprintf(fid,'double ddot_T_dCam_x;\n');

fprintf(fid,'double ddot_T_dCbm_x;\n');

fprintf(fid,'double ddot_T_du_fm;\n');

fprintf(fid,'double ddot_T_dCam_r;\n');

fprintf(fid,'double ddot_T_dCbm_r;\n');

fprintf(fid,'double ddot_T_du_em;\n');

for n = 1:1
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double ddot_T_dn1_%d_%d;\n', i, j);
        end
    end
end
for k = 1:N_c(n)
    for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
        fprintf(fid,'double ddot_T_dd1_%d;\n', j);
    end
end
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double ddot_T_dn2_%d_%d;\n', i, j);
        end
    end
end
for k = 1:N_c(n)
    for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
        fprintf(fid,'double ddot_T_dd2_%d;\n', j);
    end
end
end

for n = 1:1
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double ddot_T_ds_%d_%d;\n', i, j);
        end
    end
end
end

for k = 1:N_c(nn)
    fprintf(fid,'double ddot_T_dCx_out_1_%d;\n', k);
    fprintf(fid,'double ddot_T_dCx_out_2_%d;\n', k);
    fprintf(fid,'double ddot_T_dCr_out_1_%d;\n', k);
    fprintf(fid,'double ddot_T_dCr_out_2_%d;\n', k);
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        fprintf(fid,'double ddot_T_dC_in_%d_%d_%d;\n', n, i, k);
    end
end
end

for n = 1:nn
for k = 1:N_c(n)
    fprintf(fid,'double ddot_T_dgamma_%d_%d;\n', n, k);
end
end

if dv > 1
fprintf(fid,'double ddot_T_dgamma_feed;\n');
end
if dv > 0
fprintf(fid,'double ddot_T_dgamma_extin;\n');
fprintf(fid,'double ddot_T_dgamma_ext1;\n');
fprintf(fid,'double ddot_T_dgamma_ext2;\n');
fprintf(fid,'double ddot_T_dgamma_ext3;\n');
fprintf(fid,'double ddot_T_dgamma_rafin;\n');
fprintf(fid,'double ddot_T_dgamma_raf1;\n');
fprintf(fid,'double ddot_T_dgamma_raf2;\n');
fprintf(fid,'double ddot_T_dgamma_raf3;\n');
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double ddot_T_dC_f_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double ddot_T_dq_f_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double ddot_T_ddCdx_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'double ddot_T_ddqdx_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

for i = 1:nc
    fprintf(fid,'double ddot_T_dalpha_%d;\n', i);
end

for n = 1:nn
for s = 1:4
    fprintf(fid,'double ddot_T_dgamma_i_%d_%d;\n', n, s);
end
end

fprintf(fid,'double ddot_T_du_f;\n');

fprintf(fid,'double ddot_T_du_e;\n');

fprintf(fid,'double ddot_T_du_x;\n');

fprintf(fid,'double ddot_T_du_r;\n');

fprintf(fid,'double ddot_T_deps_ratio;\n');

for n = 1:nn
for i = 1:nc
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n))
            fprintf(fid,'double ddot_T_dC_%d_%d_%d;\n', n, i, j);
            fprintf(fid,'double ddot_T_dq_%d_%d_%d;\n', n, i, j);
        end
    end
end
end

fprintf(fid,'double ddot_T_dJ_1;\n');
fprintf(fid,'double ddot_T_dJ_2;\n');
fprintf(fid,'double ddot_T_dJ_3;\n');
fprintf(fid,'double ddot_T_dJ_4;\n');
fprintf(fid,'double ddot_T_dJ_5;\n');
fprintf(fid,'double ddot_T_dJ_6;\n');

for i = 1:nc
    fprintf(fid,'double ddot_T_dC_feed_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double ddot_T_dK_%d;\n', i);
end

fprintf(fid,'double ddot_T_dq1;\n');
fprintf(fid,'double ddot_T_dq2;\n');

for i = 1:nc
    fprintf(fid,'double ddot_T_db1_%d;\n', i);
end
for i = 1:nc
    fprintf(fid,'double ddot_T_db2_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double ddot_T_dk_h_%d;\n', i);
end

for n = 1:nn
fprintf(fid,'double ddot_T_dL_c_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double ddot_T_deps_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double ddot_T_dPeinv_%d;\n', n);
end

if dv > 1
fprintf(fid,'double ddot_T_dLfeed;\n');
end
if dv > 0
fprintf(fid,'double ddot_T_dLextin;\n');
fprintf(fid,'double ddot_T_dLext1;\n');
fprintf(fid,'double ddot_T_dLext2;\n');
fprintf(fid,'double ddot_T_dLext3;\n');
fprintf(fid,'double ddot_T_dLrafin;\n');
fprintf(fid,'double ddot_T_dLraf1;\n');
fprintf(fid,'double ddot_T_dLraf2;\n');
fprintf(fid,'double ddot_T_dLraf3;\n');
end

for s = 1:4
    fprintf(fid,'double ddot_T_du_%d;\n', s);
end

fprintf(fid,'double ddot_T_dt_s;\n');

for n = 1:nn
for i = 1:nc
    for j = 1:n_p(n)
        fprintf(fid,'double dot_dT_dC_%d_%d_%d;\n', n, i, j);
        fprintf(fid,'double dot_dT_dq_%d_%d_%d;\n', n, i, j);
    end
end
end

fprintf(fid,'double dot_dT_dJ_1;\n');
fprintf(fid,'double dot_dT_dJ_2;\n');
fprintf(fid,'double dot_dT_dJ_3;\n');
fprintf(fid,'double dot_dT_dJ_4;\n');
fprintf(fid,'double dot_dT_dJ_5;\n');
fprintf(fid,'double dot_dT_dJ_6;\n');

for i = 1:nc
    fprintf(fid,'double dot_dT_dC_feed_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double dot_dT_dK_%d;\n', i);
end

fprintf(fid,'double dot_dT_dq1;\n');
fprintf(fid,'double dot_dT_dq2;\n');

for i = 1:nc
    fprintf(fid,'double dot_dT_db1_%d;\n', i);
end
for i = 1:nc
    fprintf(fid,'double dot_dT_db2_%d;\n', i);
end

for i = 1:nc
    fprintf(fid,'double dot_dT_dk_h_%d;\n', i);
end

for n = 1:nn
fprintf(fid,'double dot_dT_dL_c_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double dot_dT_deps_%d;\n', n);
end

for n = 1:nn
fprintf(fid,'double dot_dT_dPeinv_%d;\n', n);
end

if dv > 1
fprintf(fid,'double dot_dT_dLfeed;\n');
end
if dv > 0
fprintf(fid,'double dot_dT_dLextin;\n');
fprintf(fid,'double dot_dT_dLext1;\n');
fprintf(fid,'double dot_dT_dLext2;\n');
fprintf(fid,'double dot_dT_dLext3;\n');
fprintf(fid,'double dot_dT_dLrafin;\n');
fprintf(fid,'double dot_dT_dLraf1;\n');
fprintf(fid,'double dot_dT_dLraf2;\n');
fprintf(fid,'double dot_dT_dLraf3;\n');
end

for s = 1:4
    fprintf(fid,'double dot_dT_du_%d;\n', s);
end

fprintf(fid,'double dot_dT_dt_s;\n\n');

%Building the simulation model
    idx = 0;

    for n = 1:nn
    for i = 1:nc
        for j = 1:n_p(n)
            fprintf(fid,'C_%d_%d_%d = y[%d];\n', n, i, j, idx); idx = idx + 1;
            fprintf(fid,'q_%d_%d_%d = y[%d];\n', n, i, j, idx); idx = idx + 1;
        end
    end
    end

    fprintf(fid,'J_1 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_2 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_3 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_4 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_5 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'J_6 = y[%d];\n', idx); idx = idx + 1;

    for i = 1:nc
        fprintf(fid,'C_feed_%d = y[%d];\n', i, idx); idx = idx + 1;
    end

    for i = 1:nc
        fprintf(fid,'K_%d = y[%d];\n', i, idx); idx = idx + 1;
    end

    fprintf(fid,'q1 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'q2 = y[%d];\n', idx); idx = idx + 1;

    for i = 1:nc
        fprintf(fid,'b1_%d = y[%d];\n', i, idx); idx = idx + 1;
    end
    for i = 1:nc
        fprintf(fid,'b2_%d = y[%d];\n', i, idx); idx = idx + 1;
    end

    for i = 1:nc
        fprintf(fid,'k_h_%d = y[%d];\n', i, idx); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'L_c_%d = y[%d];\n', n, idx); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'eps_%d = y[%d];\n', n, idx); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'Peinv_%d = y[%d];\n', n, idx); idx = idx + 1;
    end

    if dv > 1
    fprintf(fid,'Lfeed = y[%d];\n', idx); idx = idx + 1;
    end
    if dv > 0
    fprintf(fid,'Lextin = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lext1 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lext2 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lext3 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lrafin = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lraf1 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lraf2 = y[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'Lraf3 = y[%d];\n', idx); idx = idx + 1;
    end

    for s = 1:4
        fprintf(fid,'u_%d = y[%d];\n', s, idx); idx = idx + 1;
    end

    fprintf(fid,'t_s = y[%d];\n\n', idx);

    for n = 1:1
    fprintf(fid,'eps_ratio = (1-eps_%d)/eps_%d;\n', n, n);
    end

    fprintf(fid,'u_f = u_3 - u_2;\n');

    fprintf(fid,'u_x = u_1 - u_2;\n');

    fprintf(fid,'u_e = u_1 - u_4;\n');

    fprintf(fid,'u_r = u_3 - u_4;\n');

    for n = 1:nn
    for s = 1:4
        fprintf(fid,'gamma_i_%d_%d = u_%d*t_s/L_c_%d/eps_%d;\n', n, s, s, n, n);
    end
    end

    if dv > 1
    fprintf(fid,'gamma_feed = u_f*t_s/Lfeed/eps_%d;\n', nn);
    end
    if dv > 0
    fprintf(fid,'gamma_extin = u_1*t_s/Lextin/eps_%d;\n', nn);
    fprintf(fid,'gamma_ext1 = u_1*t_s/Lext1/eps_%d;\n', nn);
    fprintf(fid,'gamma_ext2 = u_2*t_s/Lext2/eps_%d;\n', nn);
    fprintf(fid,'gamma_ext3 = u_2*t_s/Lext3/eps_%d;\n', nn);
    fprintf(fid,'gamma_rafin = u_3*t_s/Lrafin/eps_%d;\n', nn);
    fprintf(fid,'gamma_raf1 = u_3*t_s/Lraf1/eps_%d;\n', nn);
    fprintf(fid,'gamma_raf2 = u_4*t_s/Lraf2/eps_%d;\n', nn);
    fprintf(fid,'gamma_raf3 = u_4*t_s/Lraf3/eps_%d;\n', nn);
    end

    for i = 1:nc
        fprintf(fid,'alpha_%d = k_h_%d*t_s;\n', i, i);
    end

    %Declaring the derivatives and interpolations
    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                jf = j + 1;
                fprintf(fid,'C_f_%d_%d_%d = (C_%d_%d_%d + C_%d_%d_%d)*0.5;\n', n, i, j, n, i, jf, n, i, j);
                fprintf(fid,'q_f_%d_%d_%d = (q_%d_%d_%d + q_%d_%d_%d)*0.5;\n', n, i, j, n, i, jf, n, i, j);
            end
        end
    end
    end

    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                jf = j + 1;
                fprintf(fid,'dCdx_%d_%d_%d = (C_%d_%d_%d - C_%d_%d_%d)*%d;\n', n, i, j, n, i, jf, n, i, j, ip(n));
                fprintf(fid,'dqdx_%d_%d_%d = (q_%d_%d_%d - q_%d_%d_%d)*%d;\n', n, i, j, n, i, jf, n, i, j, ip(n));
            end
        end
    end
    end

    for k = 1:N_c(1)
        for c = 0:2%(ceil(Cycles/N_c)+1)
            tkc = (c-1)*N_c(1)+k;
            if tkc >= 0
                fprintf(fid,'H_%d = (tk-%d) >= 0;\n', tkc + N_c(1), tkc);
    %             fprintf(fid,'H_%d = 1/(1+exp(-w*(t-%d)));\n', tkc + N_c(1), tkc);
            else
                fprintf(fid,'H_%d = (tk+%d) >= 0;\n', tkc + N_c(1), -tkc);
    %             fprintf(fid,'H_%d = 1/(1+exp(-w*(t+%d)));\n', tkc + N_c(1), -tkc);
            end
        end
    end

    for k = 1:N_c(1)
        fprintf(fid,'is_section1_%d = (1 - H_%d + H_%d - H_%d + H_%d);\n', k, k, k + N_c(1) - N_cx, k + N_c(1), k + N_c(1) + N_c(1) - N_cx);
        fprintf(fid,'is_section2_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cf, k + N_c(1) - N_cx, k + N_c(1) + N_c(1) - N_cf, k + N_c(1) + N_c(1) - N_cx);
        fprintf(fid,'is_section3_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cr, k + N_c(1) - N_cf, k + N_c(1) + N_c(1) - N_cr, k + N_c(1) + N_c(1) - N_cf);
        fprintf(fid,'is_section4_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k, k + N_c(1) - N_cr, k + N_c(1), k + N_c(1) + N_c(1) - N_cr);
    end
    
    for k = 1:N_c(1)
        fprintf(fid,'is_eluent_%d = (1 - H_%d + H_%d - H_%d + H_%d);\n', k, k, k + N_c(1) - 1, k + N_c(1), k + N_c(1) + N_c(1) - 1);
        fprintf(fid,'is_connection_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k, k + N_c(1) - N_cf - 1, k + N_c(1), k + N_c(1) + N_c(1) - N_cf - 1);
        fprintf(fid,'is_connection_%d = is_connection_%d + (H_%d - H_%d + H_%d - H_%d);\n', k, k, k + N_c(1) - N_cf, k + N_c(1) - 1, k + N_c(1) + N_c(1) - N_cf, k + N_c(1) + N_c(1) - 1);
        fprintf(fid,'is_feed_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cf - 1, k + N_c(1) - N_cf, k + N_c(1) + N_c(1) - N_cf - 1, k + N_c(1) + N_c(1) - N_cf);
    end
    
    %Quadrature functions terms
    for k = 1:N_c(1)
        fprintf(fid,'is_extract_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cx - 1, k + N_c(1) - N_cx, k + N_c(1) + N_c(1) - N_cx - 1, k + N_c(1) + N_c(1) - N_cx);
        fprintf(fid,'is_raffinate_%d = (H_%d - H_%d + H_%d - H_%d);\n', k, k + N_c(1) - N_cr - 1, k + N_c(1) - N_cr, k + N_c(1) + N_c(1) - N_cr - 1, k + N_c(1) + N_c(1) - N_cr);
    end

    for n = 1:nn
    for k = 1:N_c(1)
        fprintf(fid,'gamma_%d_%d = 0;\n', n, k);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_section1_%d*gamma_i_%d_1;\n', n, k, n, k, k, n);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_section2_%d*gamma_i_%d_2;\n', n, k, n, k, k, n);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_section3_%d*gamma_i_%d_3;\n', n, k, n, k, k, n);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_section4_%d*gamma_i_%d_4;\n', n, k, n, k, k, n);
        if n == nn
        if dv > 0
        if k == N_c(1)
            kf = 1;
        else
            kf = k + 1;
        end
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_extract_%d*(gamma_ext3-gamma_i_%d_1);\n', n, k, n, k, kf, n);
        fprintf(fid,'gamma_%d_%d = gamma_%d_%d + is_raffinate_%d*(gamma_raf3-gamma_i_%d_3);\n', n, k, n, k, kf, n);
        end
        end
    end
    if n == nn
    if dv > 1
    fprintf(fid,'gamma_%d_%d = gamma_feed;\n', n, N_c(nn)-6);
    end
    if dv > 0
    fprintf(fid,'gamma_%d_%d = gamma_extin;\n', n, N_c(nn)-5);
    fprintf(fid,'gamma_%d_%d = gamma_ext1;\n', n, N_c(nn)-4);
    fprintf(fid,'gamma_%d_%d = gamma_ext2;\n', n, N_c(nn)-3);
    fprintf(fid,'gamma_%d_%d = gamma_rafin;\n', n, N_c(nn)-2);
    fprintf(fid,'gamma_%d_%d = gamma_raf1;\n', n, N_c(nn)-1);
    fprintf(fid,'gamma_%d_%d = gamma_raf2;\n', n, N_c(nn));
    end
    end
    end
    
    for n = 1:nn
    nf = nn + 1 - n;
    for i = 1:nc
        for k = 1:N_c(nf)
            if n == nn
            if k == 1
                j = N_c(1)*n_p1(n);
            else
                j = (k-1)*n_p1(n);
            end
            fprintf(fid,'C_in_%d_%d_%d = 0;\n', nf, i, k);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_eluent_%d*C_%d_%d_%d*u_4/u_1;\n', nf, i, k, nf, i, k, k, n, i, j);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_connection_%d*C_%d_%d_%d;\n', nf, i, k, nf, i, k, k, n, i, j);
            if dv > 1
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_feed_%d*(u_f*C_%d_%d_%d/u_3 + u_2*C_%d_%d_%d/u_3);\n', nf, i, k, nf, i, k, k, n, i, (N_c(nn)-6)*n_p1(n), n, i, j);
            else
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_feed_%d*(u_f*C_feed_%d/u_3 + u_2*C_%d_%d_%d/u_3);\n', nf, i, k, nf, i, k, k, i, n, i, j);
            end
            else
            if dv > 1
            if k == N_c(nn)-6
            fprintf(fid,'C_in_%d_%d_%d = C_feed_%d;\n', nf, i, k, i);
            end
            end
            if dv > 0
            if k <= N_c(1)
            if k == N_c(1)
                kf = 1;
            else
                kf = k + 1;
            end
            j = k*n_p1(n);
            fprintf(fid,'C_in_%d_%d_%d = (1-is_extract_%d-is_raffinate_%d)*C_%d_%d_%d;\n', nf, i, k, kf, kf, n, i, j);
            elseif k ~= N_c(nn)-6 && k ~= N_c(nn)-5 && k ~= N_c(nn)-2
            j = (k-1)*n_p1(nf);
            fprintf(fid,'C_in_%d_%d_%d = C_%d_%d_%d;\n', nf, i, k, nf, i, j);
            elseif k ~= N_c(nn)-6
            fprintf(fid,'C_in_%d_%d_%d = 0;\n', nf, i, k);
            end
            else
            j = k*n_p1(n);
            fprintf(fid,'C_in_%d_%d_%d = C_%d_%d_%d;\n', nf, i, k, n, i, j);
            end
            end
        end
        for k = 1:N_c(nf)
            if n == nn
            if dv > 0
            if k == N_c(1)
                kf = 1;
            else
                kf = k + 1;
            end
            j = k*n_p1(nf);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_extract_%d*C_%d_%d_%d;\n', n, i, N_c(nn)-5, n, i, N_c(nn)-5, kf, nf, i, j);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_raffinate_%d*C_%d_%d_%d;\n', n, i, N_c(nn)-2, n, i, N_c(nn)-2, kf, nf, i, j);
            j = (N_c(nn)-3)*n_p1(n);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_extract_%d*C_%d_%d_%d;\n', n, i, k, n, i, k, kf, n, i, j);
            j = (N_c(nn))*n_p1(n);
            fprintf(fid,'C_in_%d_%d_%d = C_in_%d_%d_%d + is_raffinate_%d*C_%d_%d_%d;\n', n, i, k, n, i, k, kf, n, i, j);
            end
            end
        end
    end
    end
    
    %Quadrature functions terms
    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            if n == nn
            if dv > 0
            if k == N_c(nn)-3
            j = (k-1)*n_p1(n);
            fprintf(fid,'Cx_out_%d_%d = C_%d_%d_%d;\n', i, k, n, i, j);
            else
            fprintf(fid,'Cx_out_%d_%d = 0;\n', i, k);
            end
            if k == N_c(nn)
            j = (k-1)*n_p1(n);
            fprintf(fid,'Cr_out_%d_%d = C_%d_%d_%d;\n', i, k, n, i, j);
            else
            fprintf(fid,'Cr_out_%d_%d = 0;\n', i, k);
            end
            else
            if k <= N_c(1)
            if k == 1
                j = N_c(1)*n_p1(n);
            else
                j = (k-1)*n_p1(n);
            end
            fprintf(fid,'Cx_out_%d_%d = is_extract_%d*C_%d_%d_%d;\n', i, k, k, n, i, j);
            fprintf(fid,'Cr_out_%d_%d = is_raffinate_%d*C_%d_%d_%d;\n', i, k, k, n, i, j);
            end
            end
            end
        end
    end
    end

    fprintf(fid,'if(checkp || q1 > 0.0) {\n');
    for n = 1:1
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'n1_%d_%d = b1_%d*C_f_%d_%d_%d;\n', i, j, i, n, i, j);
            end
        end
    end
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'d1_%d = 1+b1_1*C_f_%d_1_%d+b1_2*C_f_%d_2_%d;\n', j, n, j, n, j);
        end
    end
    end
    fprintf(fid,'}\n');
    fprintf(fid,'if(checkp || q2 > 0.0) {\n');
    for n = 1:1
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'n2_%d_%d = b2_%d*C_f_%d_%d_%d;\n', i, j, i, n, i, j);
            end
        end
    end
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'d2_%d = 1+b2_1*C_f_%d_1_%d+b2_2*C_f_%d_2_%d;\n', j, n, j, n, j);
        end
    end
    end
    fprintf(fid,'}\n');

    %Defining the isotherms
    for n = 1:1
    for i = 1:nc
        fprintf(fid,'if(K_%d > 0.0) {\n', i);
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'s_%d_%d = alpha_%d*(K_%d*C_f_%d_%d_%d - q_f_%d_%d_%d);\n', i, j, i, i, n, i, j, n, i, j);
            end
        end
        fprintf(fid,'} else {\n');
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'s_%d_%d = -alpha_%d*q_f_%d_%d_%d;\n', i, j, i, n, i, j);
            end
        end
        fprintf(fid,'}\n');
    end
    fprintf(fid,'if(q1 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'s_%d_%d += alpha_%d*q1*n1_%d_%d/d1_%d;\n', i, j, i, i, j, j);
            end
        end
    end
    fprintf(fid,'}\n');
    fprintf(fid,'if(q2 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'s_%d_%d += alpha_%d*q2*n2_%d_%d/d2_%d;\n', i, j, i, i, j, j);
            end
        end
    end
    fprintf(fid,'}\n');
    end

    idx = 0;

    for n = 1:nn
    for i = 1:nc
        for j = 1:n_p(n)
            fprintf(fid,'dT_dC_%d_%d_%d = a[%d];\n', n, i, j, idx); idx = idx + 1;
            fprintf(fid,'dT_dq_%d_%d_%d = a[%d];\n', n, i, j, idx); idx = idx + 1;
        end
    end
    end

    fprintf(fid,'dT_dJ_1 = a[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'dT_dJ_2 = a[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'dT_dJ_3 = a[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'dT_dJ_4 = a[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'dT_dJ_5 = a[%d];\n', idx); idx = idx + 1;
    fprintf(fid,'dT_dJ_6 = a[%d];\n', idx);

    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            j = ((k-1)*n_p1(n)+1);
            fprintf(fid,'ddot_T_dfL_C_%d_%d_%d = dT_dC_%d_%d_%d*%d;\n', n, i, j, n, i, j, 2*ip(n));
            fprintf(fid,'ddot_T_dfL_q_%d_%d_%d = dT_dq_%d_%d_%d*%d;\n', n, i, j, n, i, j, 2*ip(n));
        end
    end
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                if j == ((k-1)*n_p1(n)+1)
                    fprintf(fid,'ddot_T_df_C_%d_%d_%d = dT_dC_%d_%d_%d*%d-dT_dC_%d_%d_%d*%d;\n', n, i, j, n, i, j + 1, ip(n), n, i, j, 2*ip(n));
                    fprintf(fid,'ddot_T_df_q_%d_%d_%d = dT_dq_%d_%d_%d*%d-dT_dq_%d_%d_%d*%d;\n', n, i, j, n, i, j + 1, ip(n), n, i, j, 2*ip(n));
                elseif j == (k*n_p1(n)-1)
                    fprintf(fid,'ddot_T_df_C_%d_%d_%d = dT_dC_%d_%d_%d*%d-dT_dC_%d_%d_%d*%d;\n', n, i, j, n, i, j + 1, 2*ip(n), n, i, j, ip(n));
                    fprintf(fid,'ddot_T_df_q_%d_%d_%d = dT_dq_%d_%d_%d*%d-dT_dq_%d_%d_%d*%d;\n', n, i, j, n, i, j + 1, 2*ip(n), n, i, j, ip(n));
                else              
                    fprintf(fid,'ddot_T_df_C_%d_%d_%d = (dT_dC_%d_%d_%d-dT_dC_%d_%d_%d)*%d;\n', n, i, j, n, i, j + 1, n, i, j, ip(n));
                    fprintf(fid,'ddot_T_df_q_%d_%d_%d = (dT_dq_%d_%d_%d-dT_dq_%d_%d_%d)*%d;\n', n, i, j, n, i, j + 1, n, i, j, ip(n));
                end
            end
        end
    end
    for i = 1:nc
        for k = 1:N_c(n)
            j = (k*n_p1(n));
            fprintf(fid,'ddot_T_dfR_C_%d_%d_%d = -dT_dC_%d_%d_%d*%d;\n', n, i, j, n, i, j, 2*ip(n));
            fprintf(fid,'ddot_T_dfR_q_%d_%d_%d = -dT_dq_%d_%d_%d*%d;\n', n, i, j, n, i, j, 2*ip(n));
        end
    end
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                if j == ((k-1)*n_p1(n)+1)
                    fprintf(fid,'ddot_T_ds_C_%d_%d_%d = 0.5*dT_dC_%d_%d_%d+dT_dC_%d_%d_%d;\n', n, i, j, n, i, j + 1, n, i, j);
                    fprintf(fid,'ddot_T_ds_q_%d_%d_%d = 0.5*dT_dq_%d_%d_%d+dT_dq_%d_%d_%d;\n', n, i, j, n, i, j + 1, n, i, j);
                elseif j == (k*n_p1(n)-1)
                    fprintf(fid,'ddot_T_ds_C_%d_%d_%d = dT_dC_%d_%d_%d+0.5*dT_dC_%d_%d_%d;\n', n, i, j, n, i, j + 1, n, i, j);
                    fprintf(fid,'ddot_T_ds_q_%d_%d_%d = dT_dq_%d_%d_%d+0.5*dT_dq_%d_%d_%d;\n', n, i, j, n, i, j + 1, n, i, j);
                else              
                    fprintf(fid,'ddot_T_ds_C_%d_%d_%d = 0.5*(dT_dC_%d_%d_%d+dT_dC_%d_%d_%d);\n', n, i, j, n, i, j + 1, n, i, j);
                    fprintf(fid,'ddot_T_ds_q_%d_%d_%d = 0.5*(dT_dq_%d_%d_%d+dT_dq_%d_%d_%d);\n', n, i, j, n, i, j + 1, n, i, j);
                end
            end
        end
    end
    end

    fprintf(fid,'ddot_T_dCam_x = dT_dJ_1;\n');

    fprintf(fid,'ddot_T_dCbm_x = dT_dJ_2;\n');
    
    fprintf(fid,'ddot_T_du_fm = dT_dJ_3;\n');

    fprintf(fid,'ddot_T_dCam_r = dT_dJ_4;\n');
    
    fprintf(fid,'ddot_T_dCbm_r = dT_dJ_5;\n');
    
    fprintf(fid,'ddot_T_du_em = dT_dJ_6;\n');

    for n = 1:1
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_ds_%d_%d = -eps_ratio*ddot_T_ds_C_%d_%d_%d+ddot_T_ds_q_%d_%d_%d;\n', i, j, n, i, j, n, i, j);
            end
        end
    end
    end

    for n = 1:1
    fprintf(fid,'if(q1 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dn1_%d_%d = alpha_%d*q1*ddot_T_ds_%d_%d/d1_%d;\n', i, j, i, i, j, j);
            end
        end
    end
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'ddot_T_dd1_%d = -q1*(alpha_1*n1_1_%d*ddot_T_ds_1_%d+alpha_2*n1_2_%d*ddot_T_ds_2_%d)/d1_%d/d1_%d;\n', j, j, j, j, j, j, j);
        end
    end
    fprintf(fid,'} else {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dn1_%d_%d = 0;\n', i, j);
            end
        end
    end
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'ddot_T_dd1_%d = 0;\n', j);
        end
    end
    fprintf(fid,'}\n');
    fprintf(fid,'if(q2 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dn2_%d_%d = alpha_%d*q2*ddot_T_ds_%d_%d/d2_%d;\n', i, j, i, i, j, j);
            end
        end
    end
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'ddot_T_dd2_%d = -q2*(alpha_1*n2_1_%d*ddot_T_ds_1_%d+alpha_2*n2_2_%d*ddot_T_ds_2_%d)/d2_%d/d2_%d;\n', j, j, j, j, j, j, j);
        end
    end
    fprintf(fid,'} else {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dn2_%d_%d = 0;\n', i, j);
            end
        end
    end
    for k = 1:N_c(n)
        for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
            fprintf(fid,'ddot_T_dd2_%d = 0;\n', j);
        end
    end
    fprintf(fid,'}\n');
    end

    for k = 1:N_c(nn)
        fprintf(fid,'ddot_T_dCx_out_1_%d = u_x*ddot_T_dCam_x;\n', k);
        fprintf(fid,'ddot_T_dCx_out_2_%d = u_x*ddot_T_dCbm_x;\n', k);
        fprintf(fid,'ddot_T_dCr_out_1_%d = u_r*ddot_T_dCam_r;\n', k);
        fprintf(fid,'ddot_T_dCr_out_2_%d = u_r*ddot_T_dCbm_r;\n', k);
    end

    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            j = ((k-1)*n_p1(n)+1);
            fprintf(fid,'ddot_T_dC_in_%d_%d_%d = gamma_%d_%d*ddot_T_dfL_C_%d_%d_%d;\n', n, i, k, n, k, n, i, j);
        end
    end
    end

    for n = 1:nn
    for k = 1:N_c(n)
        fprintf(fid,'ddot_T_dgamma_%d_%d = 0;\n', n, k);
    end
    for k = 1:N_c(n)
        for i = 1:nc
            j = ((k-1)*n_p1(n)+1);
            fprintf(fid,'ddot_T_dgamma_%d_%d += ddot_T_dfL_C_%d_%d_%d*C_in_%d_%d_%d;\n', n, k, n, i, j, n, i, k);
        end
    end
    fprintf(fid,'if(Peinv_%d > 0.0) {\n', n);
    for k = 1:N_c(n)
        for i = 1:nc
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dgamma_%d_%d += ddot_T_df_C_%d_%d_%d*(C_f_%d_%d_%d-dCdx_%d_%d_%d*Peinv_%d);\n', n, k, n, i, j, n, i, j, n, i, j, n);
            end
        end
    end
    fprintf(fid,'} else {\n');
    for k = 1:N_c(n)
        for i = 1:nc
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dgamma_%d_%d += ddot_T_df_C_%d_%d_%d*C_f_%d_%d_%d;\n', n, k, n, i, j, n, i, j);
            end
        end
    end
    fprintf(fid,'}\n');
    for k = 1:N_c(n)
        for i = 1:nc
            j = (k*n_p1(n));
            fprintf(fid,'ddot_T_dgamma_%d_%d += ddot_T_dfR_C_%d_%d_%d*C_%d_%d_%d;\n', n, k, n, i, j, n, i, j);
        end
    end
    end

    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dC_f_%d_%d_%d = gamma_%d_%d*ddot_T_df_C_%d_%d_%d;\n', n, i, j, n, k, n, i, j);
            end
        end
    end
    end
    for n = 1:1
    for i = 1:nc
        fprintf(fid,'if(K_%d > 0.0) {\n', i);
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dC_f_%d_%d_%d += alpha_%d*K_%d*ddot_T_ds_%d_%d;\n', n, i, j, i, i, i, j);
            end
        end
        fprintf(fid,'}\n');
    end
    fprintf(fid,'if(q1 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dC_f_%d_%d_%d += b1_%d*ddot_T_dn1_%d_%d;\n', n, i, j, i, i, j);
                fprintf(fid,'ddot_T_dC_f_%d_%d_%d += b1_%d*ddot_T_dd1_%d;\n', n, i, j, i, j);
            end
        end
    end
    fprintf(fid,'}\n');
    fprintf(fid,'if(q2 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dC_f_%d_%d_%d += b2_%d*ddot_T_dn2_%d_%d;\n', n, i, j, i, i, j);
                fprintf(fid,'ddot_T_dC_f_%d_%d_%d += b2_%d*ddot_T_dd2_%d;\n', n, i, j, i, j);
            end
        end
    end
    fprintf(fid,'}\n');
    end

    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                if n == 1
                fprintf(fid,'ddot_T_dq_f_%d_%d_%d = -alpha_%d*ddot_T_ds_%d_%d;\n', n, i, j, i, i, j);
                else
                fprintf(fid,'ddot_T_dq_f_%d_%d_%d = 0;\n', n, i, j);
                end
            end
        end
    end
    end

    for n = 1:nn
    fprintf(fid,'if(Peinv_%d > 0.0) {\n', n);
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_ddCdx_%d_%d_%d = -gamma_%d_%d*ddot_T_df_C_%d_%d_%d*Peinv_%d;\n', n, i, j, n, k, n, i, j, n);
            end
        end
    end
    fprintf(fid,'} else {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_ddCdx_%d_%d_%d = 0;\n', n, i, j);
            end
        end
    end
    fprintf(fid,'}\n');
    end

    for n = 1:nn
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_ddqdx_%d_%d_%d = 0;\n', n, i, j);
            end
        end
    end
    end

    for i = 1:nc
        fprintf(fid,'ddot_T_dalpha_%d = 0;\n', i);
    end
    for n = 1:1
    for i = 1:nc
        fprintf(fid,'if(K_%d > 0.0) {\n', i);
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dalpha_%d += ddot_T_ds_%d_%d*(K_%d*C_f_%d_%d_%d - q_f_%d_%d_%d);\n', i, i, j, i, n, i, j, n, i, j);
            end
        end
        fprintf(fid,'} else {\n');
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dalpha_%d += -ddot_T_ds_%d_%d*q_f_%d_%d_%d;\n', i, i, j, n, i, j);
            end
        end
        fprintf(fid,'}\n');
    end
    fprintf(fid,'if(q1 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dalpha_%d += ddot_T_ds_%d_%d*q1*n1_%d_%d/d1_%d;\n', i, i, j, i, j, j);
            end
        end
    end
    fprintf(fid,'}\n');
    fprintf(fid,'if(q2 > 0.0) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dalpha_%d += ddot_T_ds_%d_%d*q2*n2_%d_%d/d2_%d;\n', i, i, j, i, j, j);
            end
        end
    end
    fprintf(fid,'}\n');
    end

    for n = 1:nn
    for s = 1:4
        fprintf(fid,'ddot_T_dgamma_i_%d_%d = 0;\n', n, s);
        for k = 1:N_c(1)
            if s == 1
                fprintf(fid,'ddot_T_dgamma_i_%d_%d += is_section1_%d*ddot_T_dgamma_%d_%d;\n', n, s, k, n, k);
            elseif s == 2
                fprintf(fid,'ddot_T_dgamma_i_%d_%d += is_section2_%d*ddot_T_dgamma_%d_%d;\n', n, s, k, n, k);
            elseif s == 3
                fprintf(fid,'ddot_T_dgamma_i_%d_%d += is_section3_%d*ddot_T_dgamma_%d_%d;\n', n, s, k, n, k);
            elseif s == 4
                fprintf(fid,'ddot_T_dgamma_i_%d_%d += is_section4_%d*ddot_T_dgamma_%d_%d;\n', n, s, k, n, k);
            end
            if n == nn
            if dv > 0
            if k == N_c(1)
                kf = 1;
            else
                kf = k + 1;
            end
            if s == 1
                fprintf(fid,'ddot_T_dgamma_i_%d_%d += -is_extract_%d*ddot_T_dgamma_%d_%d;\n', n, s, kf, n, k);
            elseif s == 3
                fprintf(fid,'ddot_T_dgamma_i_%d_%d += -is_raffinate_%d*ddot_T_dgamma_%d_%d;\n', n, s, kf, n, k);
            end
            end
            end
        end
    end
    end

    for n = 1:nn
    if n == nn
    if dv > 1
    fprintf(fid,'ddot_T_dgamma_feed = ddot_T_dgamma_%d_%d;\n', n, N_c(nn)-6);
    end
    if dv > 0
    fprintf(fid,'ddot_T_dgamma_extin = ddot_T_dgamma_%d_%d;\n', n, N_c(nn)-5);
    fprintf(fid,'ddot_T_dgamma_ext1 = ddot_T_dgamma_%d_%d;\n', n, N_c(nn)-4);
    fprintf(fid,'ddot_T_dgamma_ext2 = ddot_T_dgamma_%d_%d;\n', n, N_c(nn)-3);
    fprintf(fid,'ddot_T_dgamma_ext3 = 0;\n');
    fprintf(fid,'ddot_T_dgamma_rafin = ddot_T_dgamma_%d_%d;\n', n, N_c(nn)-2);
    fprintf(fid,'ddot_T_dgamma_raf1 = ddot_T_dgamma_%d_%d;\n', n, N_c(nn)-1);
    fprintf(fid,'ddot_T_dgamma_raf2 = ddot_T_dgamma_%d_%d;\n', n, N_c(nn));
    fprintf(fid,'ddot_T_dgamma_raf3 = 0;\n');
    end
    end
    for k = 1:N_c(1)
        if n == nn
        if dv > 0
        if k == N_c(1)
            kf = 1;
        else
            kf = k + 1;
        end
        fprintf(fid,'ddot_T_dgamma_ext3 += is_extract_%d*ddot_T_dgamma_%d_%d;\n', kf, n, k);
        fprintf(fid,'ddot_T_dgamma_raf3 += is_raffinate_%d*ddot_T_dgamma_%d_%d;\n', kf, n, k);
        end
        end
    end
    end
    
    fprintf(fid,'ddot_T_du_f = ddot_T_du_fm;\n');
    for n = 1:1
    for i = 1:nc
        for k = 1:N_c(n)
            fprintf(fid,'ddot_T_du_f += is_feed_%d*ddot_T_dC_in_%d_%d_%d*C_feed_%d/u_3;\n', k, n, i, k, i);
        end
    end
    end
    if dv > 1
    fprintf(fid,'ddot_T_du_f += ddot_T_dgamma_feed*t_s/Lfeed/eps_%d;\n', nn);
    end

    fprintf(fid,'ddot_T_du_e = ddot_T_du_em;\n');
    
    fprintf(fid,'ddot_T_du_x = ddot_T_dCam_x*(');
    for k = 1:N_c(nn)%(1)
        if k < N_c(nn)%(1)
            fprintf(fid,'Cx_out_1_%d + ',k);
        else
            fprintf(fid,'Cx_out_1_%d',k);
        end
    end
    fprintf(fid,');\n');

    fprintf(fid,'ddot_T_du_x += ddot_T_dCbm_x*(');
    for k = 1:N_c(nn)%(1)
        if k < N_c(nn)%(1)
            fprintf(fid,'Cx_out_2_%d + ',k);
        else
            fprintf(fid,'Cx_out_2_%d',k);
        end
    end
    fprintf(fid,');\n');
    
    fprintf(fid,'ddot_T_du_r = ddot_T_dCam_r*(');
    for k = 1:N_c(nn)%(1)
        if k < N_c(nn)%(1)
            fprintf(fid,'Cr_out_1_%d + ',k);
        else
            fprintf(fid,'Cr_out_1_%d',k);
        end
    end
    fprintf(fid,');\n');

    fprintf(fid,'ddot_T_du_r += ddot_T_dCbm_r*(');
    for k = 1:N_c(nn)%(1)
        if k < N_c(nn)%(1)
            fprintf(fid,'Cr_out_2_%d + ',k);
        else
            fprintf(fid,'Cr_out_2_%d',k);
        end
    end
    fprintf(fid,');\n');
    
    fprintf(fid,'ddot_T_deps_ratio = 0;\n');
    fprintf(fid,'if(checkp) {\n');
    for n = 1:1
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_deps_ratio += -s_%d_%d*ddot_T_ds_C_%d_%d_%d;\n', i, j, n, i, j);
            end
        end
    end
    end
    fprintf(fid,'}\n');

    for n = 1:nn
    nf = nn + 1 - n;
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n))
                if j == ((k-1)*n_p1(n)+1)
                    fprintf(fid,'ddot_T_dC_%d_%d_%d = ddot_T_dC_f_%d_%d_%d*0.5-ddot_T_ddCdx_%d_%d_%d*%d;\n', n, i, j, n, i, j, n, i, j, ip(n));
                    fprintf(fid,'ddot_T_dq_%d_%d_%d = ddot_T_dq_f_%d_%d_%d*0.5-ddot_T_ddqdx_%d_%d_%d*%d;\n', n, i, j, n, i, j, n, i, j, ip(n));
                elseif j == (k*n_p1(n))
                    jb = j - 1;
                    fprintf(fid,'ddot_T_dC_%d_%d_%d = ddot_T_dC_f_%d_%d_%d*0.5+ddot_T_ddCdx_%d_%d_%d*%d;\n', n, i, j, n, i, jb, n, i, jb, ip(n));
                    fprintf(fid,'ddot_T_dq_%d_%d_%d = ddot_T_dq_f_%d_%d_%d*0.5+ddot_T_ddqdx_%d_%d_%d*%d;\n', n, i, j, n, i, jb, n, i, jb, ip(n));
                    if n == nn
                    if k <= N_c(1)
                    if k == N_c(1)
                        kf = 1;
                    else
                        kf = k + 1;
                    end
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += is_eluent_%d*ddot_T_dC_in_%d_%d_%d*u_4/u_1;\n', n, i, j, kf, nf, i, kf);
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += is_connection_%d*ddot_T_dC_in_%d_%d_%d;\n', n, i, j, kf, nf, i, kf);
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += is_feed_%d*u_2*ddot_T_dC_in_%d_%d_%d/u_3;\n', n, i, j, kf, nf, i, kf);
                    elseif k ~= N_c(nn)-6 && k ~= N_c(nn)-3 && k ~= N_c(nn)
                    if dv > 0
                    kf = k + 1;
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += ddot_T_dC_in_%d_%d_%d;\n', n, i, j, n, i, kf);
                    end
                    end
                    if dv > 0
                    if k == N_c(nn)-4
                    kf = k + 1;
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += ddot_T_dCx_out_%d_%d;\n', n, i, j, i, kf);
                    end
                    if k == N_c(nn)-1
                    kf = k + 1;
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += ddot_T_dCr_out_%d_%d;\n', n, i, j, i, kf);
                    end
                    else
                    if k <= N_c(1)
                    if k == N_c(1)
                        kf = 1;
                    else
                        kf = k + 1;
                    end
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += is_extract_%d*ddot_T_dCx_out_%d_%d;\n', n, i, j, kf, i, kf);
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += is_raffinate_%d*ddot_T_dCr_out_%d_%d;\n', n, i, j, kf, i, kf);
                    end
                    end
                    else
                    if dv > 0
                    if k <= N_c(1)
                    if k == N_c(1)
                        kf = 1;
                    else
                        kf = k + 1;
                    end
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += (1-is_extract_%d-is_raffinate_%d)*ddot_T_dC_in_%d_%d_%d;\n', n, i, j, kf, kf, nf, i, k);
                    end
                    else
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += ddot_T_dC_in_%d_%d_%d;\n', n, i, j, nf, i, k);
                    end
                    end
                    fprintf(fid,'ddot_T_dC_%d_%d_%d += gamma_%d_%d*ddot_T_dfR_C_%d_%d_%d;\n', n, i, j, n, k, n, i, j);
                else
                    jb = j - 1;
                    fprintf(fid,'ddot_T_dC_%d_%d_%d = (ddot_T_dC_f_%d_%d_%d+ddot_T_dC_f_%d_%d_%d)*0.5+(ddot_T_ddCdx_%d_%d_%d-ddot_T_ddCdx_%d_%d_%d)*%d;\n', n, i, j, n, i, jb, n, i, j, n, i, jb, n, i, j, ip(n));
                    fprintf(fid,'ddot_T_dq_%d_%d_%d = (ddot_T_dq_f_%d_%d_%d+ddot_T_dq_f_%d_%d_%d)*0.5+(ddot_T_ddqdx_%d_%d_%d-ddot_T_ddqdx_%d_%d_%d)*%d;\n', n, i, j, n, i, jb, n, i, j, n, i, jb, n, i, j, ip(n));
                end
            end
        end
        for k = 1:N_c(nf)
            if n == nn
            if dv > 1
            j = (N_c(nn)-6)*n_p1(n);
            fprintf(fid,'ddot_T_dC_%d_%d_%d += is_feed_%d*u_f*ddot_T_dC_in_%d_%d_%d/u_3;\n', n, i, j, k, nf, i, k);
            end
            if dv > 0
            if k == N_c(1)
                kf = 1;
            else
                kf = k + 1;
            end
            j = k*n_p1(nf);
            fprintf(fid,'ddot_T_dC_%d_%d_%d += is_extract_%d*ddot_T_dC_in_%d_%d_%d;\n', nf, i, j, kf, n, i, N_c(nn)-5);
            fprintf(fid,'ddot_T_dC_%d_%d_%d += is_raffinate_%d*ddot_T_dC_in_%d_%d_%d;\n', nf, i, j, kf, n, i, N_c(nn)-2);
            j = (N_c(nn)-3)*n_p1(n);
            fprintf(fid,'ddot_T_dC_%d_%d_%d += is_extract_%d*ddot_T_dC_in_%d_%d_%d;\n', n, i, j, kf, n, i, k);
            j = (N_c(nn))*n_p1(n);
            fprintf(fid,'ddot_T_dC_%d_%d_%d += is_raffinate_%d*ddot_T_dC_in_%d_%d_%d;\n', n, i, j, kf, n, i, k);
            end
            end
        end
    end
    end

    fprintf(fid,'ddot_T_dJ_1 = 0;\n');
    fprintf(fid,'ddot_T_dJ_2 = 0;\n');
    fprintf(fid,'ddot_T_dJ_3 = 0;\n');
    fprintf(fid,'ddot_T_dJ_4 = 0;\n');
    fprintf(fid,'ddot_T_dJ_5 = 0;\n');
    fprintf(fid,'ddot_T_dJ_6 = 0;\n');

    for i = 1:nc
        fprintf(fid,'ddot_T_dC_feed_%d = 0;\n', i);
    end
    for n = 1:1
    nf = nn + 1 - n;
    for i = 1:nc
        fprintf(fid,'if(checkp) {\n');
        if dv > 1
        k = N_c(nn)-6;
        fprintf(fid,'ddot_T_dC_feed_%d += ddot_T_dC_in_%d_%d_%d;\n', i, nf, i, k);
        else
        for k = 1:N_c(n)
            fprintf(fid,'ddot_T_dC_feed_%d += is_feed_%d*u_f*ddot_T_dC_in_%d_%d_%d/u_3;\n', i, k, n, i, k);
        end
        end
        fprintf(fid,'}\n');
    end
    end
    
    for n = 1:1
    for i = 1:nc
        fprintf(fid,'ddot_T_dK_%d = 0;\n', i);
        fprintf(fid,'if(checkp) {\n');
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dK_%d += alpha_%d*C_f_%d_%d_%d*ddot_T_ds_%d_%d;\n', i, i, n, i, j, i, j);
            end
        end
        fprintf(fid,'}\n');
    end
    end

    for n = 1:1
    fprintf(fid,'ddot_T_dq1 = 0;\n');
    fprintf(fid,'if(checkp) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dq1 += alpha_%d*ddot_T_ds_%d_%d*n1_%d_%d/d1_%d;\n', i, i, j, i, j, j);
            end
        end
    end
    fprintf(fid,'}\n');
    fprintf(fid,'ddot_T_dq2 = 0;\n');
    fprintf(fid,'if(checkp) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dq2 += alpha_%d*ddot_T_ds_%d_%d*n2_%d_%d/d2_%d;\n', i, i, j, i, j, j);
            end
        end
    end
    fprintf(fid,'}\n');
    end

    for n = 1:1
    for i = 1:nc
        fprintf(fid,'ddot_T_db1_%d = 0;\n', i);
        fprintf(fid,'if(checkp) {\n');
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_db1_%d += (ddot_T_dn1_%d_%d+ddot_T_dd1_%d)*C_f_%d_%d_%d;\n', i, i, j, j, n, i, j);
            end
        end
        fprintf(fid,'}\n');
    end
    for i = 1:nc
        fprintf(fid,'ddot_T_db2_%d = 0;\n', i);
        fprintf(fid,'if(checkp) {\n');
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_db2_%d += (ddot_T_dn2_%d_%d+ddot_T_dd2_%d)*C_f_%d_%d_%d;\n', i, i, j, j, n, i, j);
            end
        end
        fprintf(fid,'}\n');
    end
    end

    for i = 1:nc
        fprintf(fid,'ddot_T_dk_h_%d = 0;\n', i);
        fprintf(fid,'if(checkp) {\n');
        fprintf(fid,'ddot_T_dk_h_%d = ddot_T_dalpha_%d*t_s;\n', i, i);
        fprintf(fid,'}\n');
    end
    
    for n = 1:nn
    fprintf(fid,'ddot_T_dL_c_%d = 0;\n', n);
    fprintf(fid,'if(checkp) {\n');
    for s = 1:4
        fprintf(fid,'ddot_T_dL_c_%d += -u_%d*t_s/L_c_%d/L_c_%d*ddot_T_dgamma_i_%d_%d/eps_%d;\n', n, s, n, n, n, s, n);
    end
    fprintf(fid,'}\n');
    end
    
    for n = 1:nn
    fprintf(fid,'ddot_T_deps_%d = 0;\n', n);
    fprintf(fid,'if(checkp) {\n');
    for s = 1:4
        fprintf(fid,'ddot_T_deps_%d += -u_%d*t_s/L_c_%d/eps_%d/eps_%d*ddot_T_dgamma_i_%d_%d;\n', n, s, n, n, n, n, s);
    end
    if n == 1
    fprintf(fid,'ddot_T_deps_%d += -ddot_T_deps_ratio/eps_%d/eps_%d;\n', n, n, n);
    end
    if n == nn
    if dv > 1
    fprintf(fid,'ddot_T_deps_%d += -u_f*t_s/Lfeed/eps_%d/eps_%d*ddot_T_dgamma_feed;\n', n, n, n);
    end
    if dv > 0
    fprintf(fid,'ddot_T_deps_%d += -u_1*t_s/Lextin/eps_%d/eps_%d*ddot_T_dgamma_extin;\n', n, n, n);
    fprintf(fid,'ddot_T_deps_%d += -u_1*t_s/Lext1/eps_%d/eps_%d*ddot_T_dgamma_ext1;\n', n, n, n);
    fprintf(fid,'ddot_T_deps_%d += -u_2*t_s/Lext2/eps_%d/eps_%d*ddot_T_dgamma_ext2;\n', n, n, n);
    fprintf(fid,'ddot_T_deps_%d += -u_2*t_s/Lext3/eps_%d/eps_%d*ddot_T_dgamma_ext3;\n', n, n, n);
    fprintf(fid,'ddot_T_deps_%d += -u_3*t_s/Lrafin/eps_%d/eps_%d*ddot_T_dgamma_rafin;\n', n, n, n);
    fprintf(fid,'ddot_T_deps_%d += -u_3*t_s/Lraf1/eps_%d/eps_%d*ddot_T_dgamma_raf1;\n', n, n, n);
    fprintf(fid,'ddot_T_deps_%d += -u_4*t_s/Lraf2/eps_%d/eps_%d*ddot_T_dgamma_raf2;\n', n, n, n);
    fprintf(fid,'ddot_T_deps_%d += -u_4*t_s/Lraf3/eps_%d/eps_%d*ddot_T_dgamma_raf3;\n', n, n, n);
    end
    end
    fprintf(fid,'}\n');
    end

    for n = 1:nn
    fprintf(fid,'ddot_T_dPeinv_%d = 0;\n', n);
    fprintf(fid,'if(checkp) {\n');
    for i = 1:nc
        for k = 1:N_c(n)
            for j = ((k-1)*n_p1(n)+1):(k*n_p1(n)-1)
                fprintf(fid,'ddot_T_dPeinv_%d += -gamma_%d_%d*dCdx_%d_%d_%d*ddot_T_df_C_%d_%d_%d;\n', n, n, k, n, i, j, n, i, j);
            end
        end
    end
    fprintf(fid,'}\n');
    end

    if dv > 1
    fprintf(fid,'ddot_T_dLfeed = 0;\n');
    fprintf(fid,'if(checkp) {\n');
    fprintf(fid,'ddot_T_dLfeed = -u_f*t_s/Lfeed/Lfeed*ddot_T_dgamma_feed/eps_%d;\n', nn);
    fprintf(fid,'}\n');
    end
    if dv > 0
    fprintf(fid,'ddot_T_dLextin = 0;\n');
    fprintf(fid,'ddot_T_dLext1 = 0;\n');
    fprintf(fid,'ddot_T_dLext2 = 0;\n');
    fprintf(fid,'ddot_T_dLext3 = 0;\n');
    fprintf(fid,'ddot_T_dLrafin = 0;\n');
    fprintf(fid,'ddot_T_dLraf1 = 0;\n');
    fprintf(fid,'ddot_T_dLraf2 = 0;\n');
    fprintf(fid,'ddot_T_dLraf3 = 0;\n');
    fprintf(fid,'if(checkp) {\n');
    fprintf(fid,'ddot_T_dLextin = -u_1*t_s/Lextin/Lextin*ddot_T_dgamma_extin/eps_%d;\n', nn);
    fprintf(fid,'ddot_T_dLext1 = -u_1*t_s/Lext1/Lext1*ddot_T_dgamma_ext1/eps_%d;\n', nn);
    fprintf(fid,'ddot_T_dLext2 = -u_2*t_s/Lext2/Lext2*ddot_T_dgamma_ext2/eps_%d;\n', nn);
    fprintf(fid,'ddot_T_dLext3 = -u_2*t_s/Lext3/Lext3*ddot_T_dgamma_ext3/eps_%d;\n', nn);
    fprintf(fid,'ddot_T_dLrafin = -u_3*t_s/Lrafin/Lrafin*ddot_T_dgamma_rafin/eps_%d;\n', nn);
    fprintf(fid,'ddot_T_dLraf1 = -u_3*t_s/Lraf1/Lraf1*ddot_T_dgamma_raf1/eps_%d;\n', nn);
    fprintf(fid,'ddot_T_dLraf2 = -u_4*t_s/Lraf2/Lraf2*ddot_T_dgamma_raf2/eps_%d;\n', nn);
    fprintf(fid,'ddot_T_dLraf3 = -u_4*t_s/Lraf3/Lraf3*ddot_T_dgamma_raf3/eps_%d;\n', nn);
    fprintf(fid,'}\n');
    end

    for s = 1:4
        fprintf(fid,'ddot_T_du_%d = 0;\n', s);
    end
    for n = 1:nn
    nf = nn + 1 - n;
    for s = 1:4
        fprintf(fid,'ddot_T_du_%d += ddot_T_dgamma_i_%d_%d*t_s/L_c_%d/eps_%d;\n', s, n, s, n, n);
        for i = 1:nc
            for k = 1:N_c(nf)
                if n == nn
                if k == 1
                    j = N_c(1)*n_p1(n);
                else
                    j = (k-1)*n_p1(n);
                end
                if s == 1
                    fprintf(fid,'ddot_T_du_%d += -is_eluent_%d*C_%d_%d_%d*u_4/u_1/u_1*ddot_T_dC_in_%d_%d_%d;\n', s, k, n, i, j, nf, i, k);
                elseif s == 2
                    fprintf(fid,'ddot_T_du_%d += is_feed_%d*ddot_T_dC_in_%d_%d_%d*C_%d_%d_%d/u_3;\n', s, k, nf, i, k, n, i, j);
                elseif s == 3
                    fprintf(fid,'ddot_T_du_%d += -is_feed_%d*(u_f*C_feed_%d/u_3/u_3 + u_2*C_%d_%d_%d/u_3/u_3)*ddot_T_dC_in_%d_%d_%d;\n', s, k, i, n, i, j, nf, i, k);
                elseif s == 4
                    fprintf(fid,'ddot_T_du_%d += is_eluent_%d*C_%d_%d_%d*ddot_T_dC_in_%d_%d_%d/u_1;\n', s, k, n, i, j, nf, i, k);
                end
                end
            end
        end
    end
    end
    for s = 1:4
        if s == 1
            fprintf(fid,'ddot_T_du_%d += ddot_T_du_x+ddot_T_du_e;\n', s);
        elseif s == 2
            fprintf(fid,'ddot_T_du_%d += -ddot_T_du_f-ddot_T_du_x;\n', s);
        elseif s == 3
            fprintf(fid,'ddot_T_du_%d += ddot_T_du_f+ddot_T_du_r;\n', s);
        elseif s == 4
            fprintf(fid,'ddot_T_du_%d += -ddot_T_du_e-ddot_T_du_r;\n', s);
        end
    end
    if dv > 0
    for s = 1:4
        if s == 1
            fprintf(fid,'ddot_T_du_%d += ddot_T_dgamma_extin*t_s/Lextin/eps_%d+ddot_T_dgamma_ext1*t_s/Lext1/eps_%d;\n', s, nn, nn);
        elseif s == 2
            fprintf(fid,'ddot_T_du_%d += ddot_T_dgamma_ext2*t_s/Lext2/eps_%d+ddot_T_dgamma_ext3*t_s/Lext3/eps_%d;\n', s, nn, nn);
        elseif s == 3
            fprintf(fid,'ddot_T_du_%d += ddot_T_dgamma_rafin*t_s/Lrafin/eps_%d+ddot_T_dgamma_raf1*t_s/Lraf1/eps_%d;\n', s, nn, nn);
        elseif s == 4
            fprintf(fid,'ddot_T_du_%d += ddot_T_dgamma_raf2*t_s/Lraf2/eps_%d+ddot_T_dgamma_raf3*t_s/Lraf3/eps_%d;\n', s, nn, nn);
        end
    end
    end

    fprintf(fid,'ddot_T_dt_s = 0;\n');
    for i = 1:nc
        fprintf(fid,'ddot_T_dt_s += k_h_%d*ddot_T_dalpha_%d;\n', i, i);
    end
    for n = 1:nn
    for s = 1:4
        fprintf(fid,'ddot_T_dt_s += u_%d*ddot_T_dgamma_i_%d_%d/L_c_%d/eps_%d;\n', s, n, s, n, n);
    end
    end
    if dv > 1
    fprintf(fid,'ddot_T_dt_s += u_f*ddot_T_dgamma_feed/Lfeed/eps_%d;\n', nn);
    end
    if dv > 0
    for s = 1:4
        if s == 1
            fprintf(fid,'ddot_T_dt_s += u_%d*ddot_T_dgamma_extin/Lextin/eps_%d+u_%d*ddot_T_dgamma_ext1/Lext1/eps_%d;\n', s, nn, s, nn);
        elseif s == 2
            fprintf(fid,'ddot_T_dt_s += u_%d*ddot_T_dgamma_ext2/Lext2/eps_%d+u_%d*ddot_T_dgamma_ext3/Lext3/eps_%d;\n', s, nn, s, nn);
        elseif s == 3
            fprintf(fid,'ddot_T_dt_s += u_%d*ddot_T_dgamma_rafin/Lrafin/eps_%d+u_%d*ddot_T_dgamma_raf1/Lraf1/eps_%d;\n', s, nn, s, nn);
        elseif s == 4
            fprintf(fid,'ddot_T_dt_s += u_%d*ddot_T_dgamma_raf2/Lraf2/eps_%d+u_%d*ddot_T_dgamma_raf3/Lraf3/eps_%d;\n', s, nn, s, nn);
        end
    end
    end

    for n = 1:nn
    for i = 1:nc
        for j = 1:n_p(n)
            fprintf(fid,'dot_dT_dC_%d_%d_%d = -ddot_T_dC_%d_%d_%d;\n', n, i, j, n, i, j);
            fprintf(fid,'dot_dT_dq_%d_%d_%d = -ddot_T_dq_%d_%d_%d;\n', n, i, j, n, i, j);
        end
    end
    end

    fprintf(fid,'dot_dT_dJ_1 = -ddot_T_dJ_1;\n');
    fprintf(fid,'dot_dT_dJ_2 = -ddot_T_dJ_2;\n');
    fprintf(fid,'dot_dT_dJ_3 = -ddot_T_dJ_3;\n');
    fprintf(fid,'dot_dT_dJ_4 = -ddot_T_dJ_4;\n');
    fprintf(fid,'dot_dT_dJ_5 = -ddot_T_dJ_5;\n');
    fprintf(fid,'dot_dT_dJ_6 = -ddot_T_dJ_6;\n');

    for i = 1:nc
        fprintf(fid,'dot_dT_dC_feed_%d = -ddot_T_dC_feed_%d;\n', i, i);
    end

    for i = 1:nc
        fprintf(fid,'dot_dT_dK_%d = -ddot_T_dK_%d;\n', i, i);
    end

    fprintf(fid,'dot_dT_dq1 = -ddot_T_dq1;\n');
    fprintf(fid,'dot_dT_dq2 = -ddot_T_dq2;\n');

    for i = 1:nc
        fprintf(fid,'dot_dT_db1_%d = -ddot_T_db1_%d;\n', i, i);
    end
    for i = 1:nc
        fprintf(fid,'dot_dT_db2_%d = -ddot_T_db2_%d;\n', i, i);
    end

    for i = 1:nc
        fprintf(fid,'dot_dT_dk_h_%d = -ddot_T_dk_h_%d;\n', i, i);
    end
    
    for n = 1:nn
    fprintf(fid,'dot_dT_dL_c_%d = -ddot_T_dL_c_%d;\n', n, n);
    end

    for n = 1:nn
    fprintf(fid,'dot_dT_deps_%d = -ddot_T_deps_%d;\n', n, n);
    end

    for n = 1:nn
    fprintf(fid,'dot_dT_dPeinv_%d = -ddot_T_dPeinv_%d;\n', n, n);
    end

    if dv > 1
    fprintf(fid,'dot_dT_dLfeed = -ddot_T_dLfeed;\n');
    end
    if dv > 0
    fprintf(fid,'dot_dT_dLextin = -ddot_T_dLextin;\n');
    fprintf(fid,'dot_dT_dLext1 = -ddot_T_dLext1;\n');
    fprintf(fid,'dot_dT_dLext2 = -ddot_T_dLext2;\n');
    fprintf(fid,'dot_dT_dLext3 = -ddot_T_dLext3;\n');
    fprintf(fid,'dot_dT_dLrafin = -ddot_T_dLrafin;\n');
    fprintf(fid,'dot_dT_dLraf1 = -ddot_T_dLraf1;\n');
    fprintf(fid,'dot_dT_dLraf2 = -ddot_T_dLraf2;\n');
    fprintf(fid,'dot_dT_dLraf3 = -ddot_T_dLraf3;\n');
    end

    for s = 1:4
        fprintf(fid,'dot_dT_du_%d = -ddot_T_du_%d;\n', s, s);
    end

    fprintf(fid,'dot_dT_dt_s = -ddot_T_dt_s;\n\n');

    idx = 0;

    % fprintf(fid,'dadt = zeros(size(a));\n');

    for n = 1:nn
    for i = 1:nc
        for j = 1:n_p(n)
            fprintf(fid,'dadt[%d] = dot_dT_dC_%d_%d_%d;\n', idx, n, i, j); idx = idx + 1;
            fprintf(fid,'dadt[%d] = dot_dT_dq_%d_%d_%d;\n', idx, n, i, j); idx = idx + 1;
        end
    end
    end

    fprintf(fid,'dadt[%d] = dot_dT_dJ_1;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dJ_2;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dJ_3;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dJ_4;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dJ_5;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dJ_6;\n', idx); idx = idx + 1;

    for i = 1:nc
        fprintf(fid,'dadt[%d] = dot_dT_dC_feed_%d;\n', idx, i); idx = idx + 1;
    end

    for i = 1:nc
        fprintf(fid,'dadt[%d] = dot_dT_dK_%d;\n', idx, i); idx = idx + 1;
    end

    fprintf(fid,'dadt[%d] = dot_dT_dq1;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dq2;\n', idx); idx = idx + 1;

    for i = 1:nc
        fprintf(fid,'dadt[%d] = dot_dT_db1_%d;\n', idx, i); idx = idx + 1;
    end
    for i = 1:nc
        fprintf(fid,'dadt[%d] = dot_dT_db2_%d;\n', idx, i); idx = idx + 1;
    end

    for i = 1:nc
        fprintf(fid,'dadt[%d] = dot_dT_dk_h_%d;\n', idx, i); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'dadt[%d] = dot_dT_dL_c_%d;\n', idx, n); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'dadt[%d] = dot_dT_deps_%d;\n', idx, n); idx = idx + 1;
    end

    for n = 1:nn
    fprintf(fid,'dadt[%d] = dot_dT_dPeinv_%d;\n', idx, n); idx = idx + 1;
    end

    if dv > 1
    fprintf(fid,'dadt[%d] = dot_dT_dLfeed;\n', idx); idx = idx + 1;
    end
    if dv > 0
    fprintf(fid,'dadt[%d] = dot_dT_dLextin;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dLext1;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dLext2;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dLext3;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dLrafin;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dLraf1;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dLraf2;\n', idx); idx = idx + 1;
    fprintf(fid,'dadt[%d] = dot_dT_dLraf3;\n', idx); idx = idx + 1;
    end

    for s = 1:4
        fprintf(fid,'dadt[%d] = dot_dT_du_%d;\n', idx, s); idx = idx + 1;
    end

    fprintf(fid,'dadt[%d] = dot_dT_dt_s;\n\n', idx);

% fprintf(fid,'end\n');

fprintf(fid,'}\n');

fclose(fid);

% dt = 0.001;
% reps = ceil(Cycles/(tf-t0));
% nt = (tf-t0)/dt;
% eval(['codegen matlab_SMB_odefun_y -args {0,zeros(',num2str(ny),',1),0}']);
% eval(['codegen matlab_SMB_odefun_a -args {0,zeros(',num2str(ny),',1),zeros(1,',num2str(reps*nt+1),'),zeros(',num2str(ny),',',num2str(reps*nt+1),'),0}']);
% eval(['codegen -c -config:lib matlab_SMB_fun_y -args {0,zeros(',num2str(ny),',1)}']);
% eval(['codegen -c -config:lib matlab_SMB_fun_a -args {0,zeros(',num2str(ny),',1),zeros(1,',num2str(ny),')}']);
% mex -output matlab_SMB_odefun_y_mex ode_mex.c rp45.c ntrp45.c matlab_SMB_eventfun.c matlab_SMB_odefun_y.c matlab_SMB_fun_y.c
% mex -output matlab_SMB_odefun_a_mex ode_mex.c rp45.c ntrp45.c matlab_SMB_eventfun.c matlab_SMB_odefun_a.c matlab_SMB_fun_a.c binterp.c
mex -output matlab_SMB_odefun_2_y_mex ode_y_mex.c matlab_SMB_eventfun.c matlab_SMB_odefun_2_y.c matlab_SMB_fun_y.c
mex -output matlab_SMB_odefun_2_a_mex ode_a_mex.c matlab_SMB_eventfun.c matlab_SMB_odefun_2_a.c matlab_SMB_fun_a.c
mex -output matlab_SMB_fun_y_mex matlab_SMB_fun_y_mex.c matlab_SMB_fun_y.c
mex -output matlab_SMB_fun_a_mex matlab_SMB_fun_a_mex.c matlab_SMB_fun_a.c
end
