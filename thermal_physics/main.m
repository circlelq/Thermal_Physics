global N temp L sys time_step n_Itr direc_vec storage_Ek storage_Potential storage_a storage_temp chem_p boltz count nsample

%�����趨
L = 5;  %�߳�
N = 27; %������
temp = 5; %�¶�
time_step = 0.5e-3;
n_Itr = 9000; %��������
nsample = 5; %�������
count = 0; %��¼��������
boltz = 0;
%������������
obj_particles(N,1) = Particle;
sys = obj_particles; 


%���ݴ���
%�±�1������������ �±�2�� ����ID
storage_Ek = zeros(n_Itr,1);
storage_Potential = zeros(n_Itr,1);
storage_a = zeros(n_Itr,N,3);
storage_temp = zeros(n_Itr,1);
chem_p = zeros(n_Itr,1);

%����������������������
%�±�1�� ������ �±�2�� ����ID
direc_vec = zeros(3,27);

%���㷨
create_particles();
verlet_loop();
post_process();



function create_particles()
global N temp L sys

%��ʼ��λ��
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

%��ʼ���ٶ�
for i =1:N
    % ��ʼ�ٶ���(-0.5,0.5)�Ͼ��ȷֲ�
    rand_vel = zeros(3,1);
    rand_vel(1,1) = rand() - 0.5 ;
    rand_vel(2,1) = rand() - 0.5 ;
    rand_vel(3,1) = rand() - 0.5 ;
    sys(i,1).v = rand_vel;
end

%������ʼ������ʹ�ó�ʼ�����غ�
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
%��ѭ����verlet�㷨��������ϵͳ�ݻ�
global L n_Itr direc_vec nsample count

%���ɷ�������
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

for i = 1 : 100 %ǰ100�����б궨
    velocity_calibration(); %�Ƚ����ٶȱ궨
    verlet(i);
end

%����1000��
for i = 101:1000
    verlet(i);
end

%ƽ���ʼ�������㻯ѧ��
for i = 1001:n_Itr
    verlet(i);
    if mod(i,nsample)==0
        count = count + 1;
        Widom(i);
    end
end

end

function velocity_calibration()
global N temp sys

%�ٶȱ궨
instant_temp = 0; %˲ʱ�¶�
for i = 1:N
    instant_temp = instant_temp + sys(i,1).Ek*2;
end
instant_temp = instant_temp/3/N;
factor = sqrt(temp / instant_temp);
for i = 1:N
    sys(i,1).v = sys(i,1).v .* factor;
end

end



function post_process()
global time_step n_Itr storage_Ek storage_Potential storage_temp nsample chem_p count

figure(1)
tspan = 0:time_step:time_step*(n_Itr-1);
plot(tspan,storage_Ek);
hold on
plot(tspan,storage_Potential);
plot(tspan,storage_Potential+storage_Ek);
legend('Ek','Potential','total');
xlabel('time (s)');
ylabel('Energy');
title('Energy - Time Relation');
hold off

figure(2)
tspan1 = 0:time_step:time_step*(n_Itr-1);
plot(tspan1,storage_temp);
legend('Temperature');
xlabel('time (s)');
ylabel('Temperature');
title('Temperature - Time Relation');

figure(3)
tspan2 = (1000+99*nsample)*time_step:nsample*time_step:time_step*(n_Itr-1);
plot(tspan2,chem_p(100:count,1));
legend('Chemical Potential');
xlabel('time (s)');
ylabel('Chemical Potential');
title('Chemical Potential - Time Relation');
end
