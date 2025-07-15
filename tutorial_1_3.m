Ap = [1 1;0 1];
Bp = [0.5; 1];
Cp = [1 0];
Dp = 0;
Np = 70;
Nc = 4;
R_wei_var = 1;
R_wei = R_wei_var*eye(Nc);

[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e] =mpcgain(Ap,Bp,Cp,Nc,Np);

[n, n_in] = size(B_e);
xm=[0;0];
Xf = zeros(n,1);
N_sim = 100;
r=ones(N_sim,1);

u=0;
y=0;

for kk=1:N_sim

    DeltaU= inv(Phi_Phi+R_wei)*(Phi_R.*r(kk)-Phi_F*Xf);
    deltau=DeltaU(1,1); 
    u=u+deltau; %u(k)=u(k-1) + deltau
    u1(kk) =u;
    y1(kk)=y;

    xm_old = xm;
    xm = Ap*xm+Bp*u;
    y=Cp*xm;
    Xf=[xm-xm_old;y];

end

%%%%%  Closed Loop Control System %%%%%
Ky_I = zeros(1,Nc);
Ky_I(1) = 1;
Ky_MPC = zeros(1,Nc);
Ky_MPC(1) = 1;

Ky = Ky_I*inv(Phi_Phi+R_wei)*Phi_R;
Kmpc = Ky_MPC*inv(Phi_Phi+R_wei)*Phi_F;

pole_mpc = eig(A_e-B_e*Kmpc);

H = Phi_Phi+R_wei;
con_num = cond(H);

k=0:(N_sim-1);
figure(1)
subplot(2,1,1);
plot(k,y1);
grid on;
xlabel('Sampling Instant')
legend('Output')

subplot(2,1,2)
plot(k,u1)
grid on;
xlabel('Sampling Instant')
legend('Control')