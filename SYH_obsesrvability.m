% 수치값 지정
gamma   = 0;
Vx      = 0.01;
Vy      = 0.0001;
delta_f = 0.1;
delta_r = 0.05;
L       = 1.0;
M       = 50;
I       = 10;

omega_f = 10;
omega_r = 10;
r       = 0.3;
Cx      = 1000;
Cy      = 1500;
tau_x   = 0.3;
tau_y   = 0.4;

% A matrix 수치 계산
A = [
  0, gamma, Vy, cos(delta_f)/M, -sin(delta_f)/M, cos(delta_r)/M, -sin(delta_r)/M;
 -gamma, 0, -Vx, sin(delta_f)/M, cos(delta_f)/M, sin(delta_r)/M, cos(delta_r)/M;
  0, 0, 0, (sin(delta_f)*L)/(2*I), (cos(delta_f)*L)/(2*I), -(sin(delta_r)*L)/(2*I), -(cos(delta_r)*L)/(2*I);
 -Cx/(tau_x*omega_f*r), 0, 0, -1/tau_x, 0, 0, 0;
 (Cy*(Vy + L/2*gamma))/(tau_y*Vx^2), -Cy/(tau_y*Vx), -Cy*L/(2*tau_y*Vx), 0, -1/tau_y, 0, 0;
 -Cx/(tau_x*omega_r*r), 0, 0, 0, 0, -1/tau_x, 0;
 (Cy*(Vy - L/2*gamma))/(tau_y*Vx^2), -Cy/(tau_y*Vx), Cy*L/(2*tau_y*Vx), 0, 0, 0, -1/tau_y
];

% 결과 출력
disp('Evaluated A matrix:');
disp(A);

Num_s = 7;
C = [0 0 1 1 1 1 1];

Ob = C;

for i = 1:Num_s-1
    Ob = [Ob; C * A^i];
end

rank_Ob = rank(Ob);



