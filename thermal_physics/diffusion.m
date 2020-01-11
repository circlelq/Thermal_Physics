function [diff_vel,diff_pos] = diffusion()
%id_Itr: ѭ���Ĳ���
global N time_step storage_velocity storage_position n_Itr

% �����ٶȹ���
t_up = 3e3;
M = 1e3;
t0 = n_Itr - 5e3;
temp1 = 0;
for ti = 0:t_up  % ��Ӧt
    for k = 1:M   % ��Ӧ t0
        for i = 1:N  % ��Ӧÿ������
            temp1 = temp1 + storage_velocity(i,t0+k,1)*storage_velocity(i,t0+k+ti,1)...
                + storage_velocity(i,t0+k,2)*storage_velocity(i,t0+k+ti,2) + ...
                storage_velocity(i,t0+k,3)*storage_velocity(i,t0+k+ti,3);
        end
    end
end
diff_vel = temp1 * time_step / (3*M*N);

% ����λ�ù���
temp2 = 0;
for t = 1001:1100
    for k = 1:M
        for i = 1:N
            temp2 = temp2 + ...
                dot((storage_position(i,2*t+1+t0+k,:)-storage_position(i,t0+k,:))...
                ,(storage_position(i,2*t+1+t0+k,:)-storage_position(i,t0+k,:))) - ...
                dot((storage_position(i,2*t+t0+k,:)-storage_position(i,t0+k,:))...
                ,(storage_position(i,2*t+t0+k,:)-storage_position(i,t0+k,:)));
        end
    end
end
diff_pos = temp2 / (6*M*N*time_step*100);
end

