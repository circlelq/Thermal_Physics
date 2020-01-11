global N temp L sys  n_Itr direc_vec storage_Potential potential potential_old chem_p boltz count

%�����趨
L = 5;  %�߳�
N = 27; %������
temp = 5; %�¶�
n_Itr = 1e4; %����
nsample = 5; %�������
count = 0; %��¼��������
boltz = 0;
%������������
obj_particles(N,1) = Particle;
sys = obj_particles;

%���ݴ���
%�±�1����������
storage_Potential = zeros(n_Itr,1);
chem_p = zeros(n_Itr,1);
%����������������������
%�±�1�� ������ �±�2�� ����ID
direc_vec = zeros(3,27);


%���㷨

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
create_particles();
calculate_potential();
potential_old = potential;
% ����1000��
for i = 1:1000
    move_a_particle(i);
end
% ƽ���ʼ����
for i = 1001:n_Itr
    move_a_particle(i);
    if mod(i,nsample)==0
        count = count + 1;
        widom();
    end        
end
tspan = 1000+100*nsample:nsample:n_Itr;
plot(tspan,chem_p(100:count,1));
xlabel('step');
ylabel('chemical potential');
legend('Chemical Potential');
title('Chemical Potential - Step Relation');




function create_particles()
global L sys

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

end


function calculate_potential()

global N L sys direc_vec potential

potential = 0;

    for i = 1:N
        for j = i+1:N
            for k = 1 : 27  %����������
                relative_R = sys(i,1).pos - sys(j,1).pos + direc_vec(:,k);
                R = sum(relative_R.*relative_R); %�������
                term1 = R^(-6);
                term2 = R^(-3);
                
                if R < 0.25*L*L %С�ڽضϾ���
                    %��������
                    potential = potential + 4*(term1 - term2);
                    break
                else
                    continue
                end
            end
        end
    end
end


function move_a_particle(id_Itr)

global N temp L sys storage_Potential potential potential_old

delta = 1;  %ÿ�������ϵ����ƶ���������

index = unidrnd(N);  %���ѡһ�������ƶ�
position_old = sys(index,1).pos;  % �洢�ƶ�ǰ��λ��
sys(index,1).pos = sys(index,1).pos + (rand(3,1)-0.5)*delta;

%���ڱ߽�����
if sys(index,1).pos(1,1)>L || sys(index,1).pos(1,1)<0
    sys(index,1).pos(1,1) = mod(sys(index,1).pos(1,1),L);
end
if sys(index,1).pos(2,1)>L || sys(index,1).pos(2,1)<0
    sys(index,1).pos(2,1) = mod(sys(index,1).pos(2,1),L);
end
if sys(index,1).pos(3,1)>L || sys(index,1).pos(3,1)<0
    sys(index,1).pos(3,1) = mod(sys(index,1).pos(3,1),L);
end

calculate_potential();
dU = potential - potential_old;

% �����������˻�ԭλ
if rand > exp(-dU/temp)
    sys(index,1).pos = position_old;
    storage_Potential(id_Itr,1) = potential_old;
else
    storage_Potential(id_Itr,1) = potential;
    potential_old = potential;
end

end


function widom()
global N temp L sys direc_vec boltz chem_p count
new_particle = Particle;
new_particle.pos = (rand(3,1)-0.5)*L;
new_potential = 0;
for i = 1:N
    for k = 1:27
        relative_R = new_particle.pos - sys(i,1).pos + direc_vec(:,k);
        R = sum(relative_R.*relative_R); %�������
        term1 = R^(-6);
        term2 = R^(-3);
                
        if R < 0.25*L*L %С�ڽضϾ���
            %��������
            new_potential = new_potential + 4*(term1 - term2);
            break
        else
            continue
        end
    end
end
boltz = boltz + exp(-new_potential/temp);
chem_p(count,1) = -temp*log(boltz/count);
end


