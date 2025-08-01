Jm = 0.01;
Bm = 0.001;
Jl = 1;
Bl = 2;
ks = 7000;
N = 100;

%SISO System
Ac = [0 1 0 0;-ks/(N^2*Jm) -Bm/Jm ks/(Jm*N) 0; 0 0 0 1;ks/(Jl*N) 0 -ks/Jl -Bl/Jl;];
Bc = [0; 0; 0; 1;];
Cc = [1 0 0 0];
Dc = 0;
dt = 0.01;
Nc = 4;
Np= 1000;
R_wei_var = 1;
R_wei = R_wei_var*eye(Nc);

[Ap,Bp,Cp,Dp]=c2dm(Ac,Bc,Cc,Dc,dt);
[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e] =mpcgain(Ap,Bp,Cp,Nc,Np);

%%%%% Observable Check %%%%%
Num_s = 4;
Ob = Cp;
for i = 1:Num_s-1
    Ob = [Ob; Cp * Ap^i];
end
rank_Ob = rank(Ob);

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

%%%%% Initia  Value %%%%%
[n,n_in] = size(B_e);
xm = [0.05;-0.01;0;0];
N_sim = 10000;
r = 0.1*ones(N_sim,1);
Xf = zeros(n,1);

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

k=0:(N_sim-1);
figure(1)
subplot(2,1,1);
plot(k,y1);
grid on;
xlabel('Sampling Instant')
ylabel('Position');
legend('Output')

subplot(2,1,2)
plot(k,u1)
grid on;
xlabel('Sampling Instant')
ylabel('Control Input');
legend('Control')

