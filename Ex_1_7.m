
omega=2;

Ac = [0 1; -omega^2 0];
Bc = [0;1];
Cc=[0 1];
Dc = 0;
dt = 0.1;
N_sim = 100;

[Ad,Bd,Cd,Dd] = c2dm(Ac,Bc,Cc,Dc,dt);

u=0;
y=0;
xm=[1;0];
xhat = [0.3;0];

for kk=1:N_sim

    u1(kk) =u;
    y1(kk)=y;
    x1(kk)=xm(1);
    x2(kk)=xm(2);
    x1_hat(kk) = xhat(1);
    x2_hat(kk) = xhat(2);

%     xm_old = xm;
    xm = Ad*xm+Bd*u;
    y=Cd*xm;

%     xhat = Ad*xhat+Bd*u;
%     y_hat=Cd*xhat;
    xhat = Ad*xhat+Bd*u + [-1.6284;1.6601]*(Cd*xm-Cd*xhat);
    y_hat=Cd*xhat;
%     Xf=[xm-xm_old;y];

end

k=0:(N_sim-1);
figure(1)
subplot(2,1,1);
plot(k, x1,'k');
hold on
plot(k, x1_hat,'r');
grid on;
xlabel('Sampling Instant')
ylabel('\theta')
legend('true','observer');

subplot(2,1,2)
plot(k, x2,'k');
hold on
plot(k, x2_hat,'r');
grid on;
xlabel('Sampling Instant')
ylabel('\dot{\theta}')
legend('true','observer');



