%% init
clear; close all; clc;
restoredefaultpath;
gait_path = 'gait';
addpath(gait_path);
casia_gait = load('gait/casia_gait_0.mat');

%% param init
transmission_range = 50;
motor_type = 1;
figure();
for t=1:3
    motor_type = t;
switch motor_type
    case 1 %TQ 8526
        Kt = 0.127;
        Km = 0.254;
        Im = 0.000115;
        bm = 0;
        Tm_max = 8.3;
        joint_index = [9,10,16,17];
        u_index = [3,4,8,9];
    case 2 %7615
        Kt = 0.131;
        Tm_max = 2.86;
        joint_index = [7,8,14,15];
        u_index = [1,2,6,7];
    case 3 %6013
        Kt = 0.11;
        Tm_max = 1.37;
        joint_index = [13,20];
        u_index = [5,10];
end 

%% compute load rms
u_all = [casia_gait.gait(1).inputs.u, casia_gait.gait(3).inputs.u];
dx_all = [casia_gait.gait(1).states.dx, casia_gait.gait(3).states.dx];
ddx_all = [casia_gait.gait(1).states.ddx, casia_gait.gait(3).states.ddx];

u = u_all(u_index, :);
dx = dx_all(joint_index, :);
ddx = ddx_all(joint_index, :);

[TL_max, TL_max_dx_index] = max(abs(u),[],2);
TL_max_dx = [];
for i = 1:length(joint_index)
    TL_max_dx = [TL_max_dx, dx(i, TL_max_dx_index(i))];
end
TL_max_dx = TL_max_dx';

TL_rms = rms(u,2);
ddx_rms = rms(ddx,2);

%% compute total E
C1 = sum(u.*dx, 2);
C2 = sum(dx.*ddx, 2);
C3 = sum(dx.*dx, 2);
C4 = sum(ddx.*ddx, 2);
C5 = sum(u.*u, 2);
C6 = sum(u.*ddx, 2);
C7 = sum(u.*dx, 2);
E_tot_all = [];
for belta = 1:1:transmission_range
    E_tot = C1 + Im*(1+2*bm/Km^2)*belta^2*C2 +...
        bm*(1+bm/Km^2)*belta^2*C3 + ...
        Im^2/Km^2*belta^2*C4 + 1/(Km^2*belta^2)*C5 +...
        2*Im/Km^2*C6 + 2*bm/Km^2*C7;    
    E_tot_all = [E_tot_all, E_tot];
end

%% compute bandwidth
w_all = [];
for belta = 1:1:transmission_range
    w = sqrt((Tm_max/belta-TL_max/belta^2)./(TL_max_dx.*TL_max_dx)-bm^2)/Im;
    w_all = [w_all, w];
end

plot(w_all(1,:),E_tot_all(1,:),'LineWidth',2)

hold on

end
legend('8526','7615','6013')
xlabel('band width')
ylabel('energy')

xlim([0 1.2e+05])





