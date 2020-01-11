global N temp L sys time_step n_Itr direc_vec storage_Ek storage_Potential ...
    storage_a storage_temp storage_velocity storage_position...
    chem_p boltz count nsample storage_pressure

%²ÎÊýÉè¶¨
L = 5;  %±ß³¤
N = 27; %Á£×ÓÊý
temp = 5; %ÎÂ¶È
time_step = 0.5e-3;
n_Itr = 5e4; %µü´ú²½Êý
nsample = 5; %²ÉÑù¼ä¸ô
count = 0; %¼ÇÂ¼²ÉÑù´ÎÊý
boltz = 0;
%´´½¨¶ÔÏóÊý×é
obj_particles(N,1) = Particle;
sys = obj_particles; 


%Êý¾Ý´¢´æ
%ÏÂ±ê1£ºµü´ú²½Êý£» ÏÂ±ê2£º Á£×ÓID
storage_Ek = zeros(n_Itr,1);
storage_Potential = zeros(n_Itr,1);
storage_a = zeros(n_Itr,N,3);
storage_temp = zeros(n_Itr,1);
storage_pressure = zeros(n_Itr,1);
storage_velocity = zeros(N,n_Itr,3);
storage_position = zeros(N,n_Itr,3);
diffCoef_vel = 0; % ËÙ¶È¹ØÁª
diffCoef_pos = 0; % Î»ÖÃ¹ØÁª
chem_p = zeros(n_Itr,1);


%·½ÏòÏòÁ¿£¬ÓÃÓÚËÑË÷Á£×Ó
%ÏÂ±ê1£º ·ÖÁ¿£» ÏÂ±ê2£º ·½ÏòID
direc_vec = zeros(3,27);

%Ö÷Ëã·¨
create_particles();
verlet_loop();

% Ñ­»·Íâ¸üÐÂÑ¹Ç¿

storage_pressure = storage_pressure /3 /L^3 + N * storage_temp / L^3 ...
    + 16*pi/3 * (N/L^3)^2 * (2/3/(0.5*L)^9 - 1/(2*L)^3); % ÐÞÕýÏî

% ¼ÆËãÀ©É¢ÏµÊý
[diffCoef_vel,diffCoef_pos] = diffusion();

%%
% »­Í¼
post_process();

% ±£´æÎÄ¼þ
save data_test5

function create_particles()
global N temp L sys

%³õÊ¼»¯Î»ÖÃ
position = zeros(3,1);
for i = 1:3
    for j = 1:3
        for k = 1:3
            position(1,1) = (i-1)*L/3 + L/6;
            position(2,1) = (j-1)*L/3 + L/6;
            position(3,1) = (k-1)*L/3 + L/6;
            index = (i-1)*9 + (j-1)*3 +k;
            sys(index,1).pos = position;
        end
    end
end

%³õÊ¼»¯ËÙ¶È
for i =1:N
    % ³õÊ¼ËÙ¶ÈÔÚ(-0.5,0.5)ÉÏ¾ùÔÈ·Ö²¼
    rand_vel = zeros(3,1);
    rand_vel(1,1) = rand() - 0.5 ;
    rand_vel(2,1) = rand() - 0.5 ;
    rand_vel(3,1) = rand() - 0.5 ;
    sys(i,1).v = rand_vel;
end

%µ÷Õû³õÊ¼¶¯Á¿£¬Ê¹µÃ³õÊ¼¶¯Á¿ÊØºã
res_moment = zeros(3,1);
for i = 1:N
    res_moment = res_moment + sys(i,1).v;
end
res_moment = res_moment/N;
for i = 1:N
    sys(i,1).v = sys(i,1).v - res_moment;
end

end


function verlet_loop()
%Ö÷Ñ­»·£¬verletËã·¨¼ÆËãÁ£×ÓÏµÍ³ÑÝ»¯
global L n_Itr direc_vec nsample count Ut storage_Potential

%Éú³É·½ÏòÏòÁ¿
for i = 1:3
    for j = 1:3
        for k = 1:3
            index = (i-1)*9 + (j-1)*3 + k;
            direc_vec(1,index) = (i-2)*L;
            direc_vec(2,index) = (j-2)*L;
            direc_vec(3,index) = (k-2)*L;
        end
    end
end

%½Ø¶ÏÊÆÄÜ
Ut = 1 / 4 * (1/(0.5*L)^12 - (0.5*L)^6);

for i = 1 : 1000 %Ç°100²½½øÐÐ±ê¶¨
    velocity_calibration(); %ÏÈ½øÐÐËÙ¶È±ê¶¨
    verlet(i);
end

%ÏÈËã1000²½
for i = 1001:2000
    verlet(i);
    pressure(i);

end

%Æ½ºâºó¿ªÊ¼²ÉÑù¼ÆËã»¯Ñ§ÊÆ
for i = 2001:n_Itr
    verlet(i);
    if mod(i,nsample)==0
        count = count + 1;
        Widom(i);
    end
    pressure(i);
end

storage_Potential = storage_Potential - Ut;

end

function velocity_calibration()
global N temp sys

%ËÙ¶È±ê¶¨
instant_temp = 0; %Ë²Ê±ÎÂ¶È
for i = 1:N
    instant_temp = instant_temp + sys(i,1).Ek*2;
end
instant_temp = instant_temp/3/N;
factor = sqrt((temp - 0.03) / instant_temp);
for i = 1:N
    sys(i,1).v = sys(i,1).v .* factor;
end

end



function post_process()
global time_step n_Itr storage_Ek storage_Potential storage_temp nsample chem_p count storage_pressure

subplot(2,2,1)
tspan = 0:time_step:time_step*(n_Itr-1);
plot(tspan,storage_Ek);
hold on
plot(tspan,storage_Potential);
plot(tspan,storage_Potential+storage_Ek);
legend('Ek','Potential','total');
xlabel('time');
ylabel('Energy');
title('Energy - Time Relation');
hold off

subplot(2,2,2)
tspan1 = 0:time_step:time_step*(n_Itr-1);
plot(tspan1,storage_temp);
legend('Temperature');
xlabel('time');
ylabel('Temperature');
title('Temperature - Time Relation');

subplot(2,2,3)
tspan2 = (2000+99*nsample)*time_step:nsample*time_step:time_step*(n_Itr-1);
plot(tspan2,chem_p(100:count,1));
legend('Chemical Potential');
xlabel('time');
ylabel('Chemical Potential');
title('Chemical Potential - Time Relation');

subplot(2,2,4)
tspan1 = 0:time_step:time_step*(n_Itr-1);
plot(tspan1,storage_pressure');
legend('Pressure');
xlabel('time');
ylabel('Pressure');
title('Pressure - Time Relation');
end
