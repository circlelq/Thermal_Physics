global N temp L sys  n_Itr direc_vec storage_Potential potential potential_old chem_p boltz count

%参数设定
L = 5;  %边长
N = 27; %粒子数
temp = 5; %温度
n_Itr = 1e4; %步数
nsample = 5; %采样间隔
count = 0; %记录采样次数
boltz = 0;
%创建对象数组
obj_particles(N,1) = Particle;
sys = obj_particles;

%数据储存
%下标1：迭代步数
storage_Potential = zeros(n_Itr,1);
chem_p = zeros(n_Itr,1);
%方向向量，用于搜索粒子
%下标1： 分量； 下标2： 方向ID
direc_vec = zeros(3,27);


%主算法

%生成方向向量
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
% 先算1000步
for i = 1:1000
    move_a_particle(i);
end
% 平衡后开始采样
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

%初始化位置
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
            for k = 1 : 27  %搜索像粒子
                relative_R = sys(i,1).pos - sys(j,1).pos + direc_vec(:,k);
                R = sum(relative_R.*relative_R); %计算距离
                term1 = R^(-6);
                term2 = R^(-3);
                
                if R < 0.25*L*L %小于截断距离
                    %计算势能
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

delta = 1;  %每个方向上单次移动的最大距离

index = unidrnd(N);  %随机选一个粒子移动
position_old = sys(index,1).pos;  % 存储移动前的位置
sys(index,1).pos = sys(index,1).pos + (rand(3,1)-0.5)*delta;

%周期边界条件
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

% 若不接受则退回原位
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
        R = sum(relative_R.*relative_R); %计算距离
        term1 = R^(-6);
        term2 = R^(-3);
                
        if R < 0.25*L*L %小于截断距离
            %计算势能
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


