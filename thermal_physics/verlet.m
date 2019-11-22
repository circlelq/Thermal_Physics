function verlet(id_Itr)
%������verlet�㷨
%id_Itr: ��ǰ��������

global N L sys direc_vec time_step storage_a storage_Ek storage_Potential storage_temp n_Itr

%�����ٽ����Ӽ�������
if id_Itr == 1
    for i = 1:N
        for j = i+1:N
            for k = 1 : N  %����������
                relative_R = sys(i,1).pos - sys(j,1).pos + direc_vec(:,k);
                R = sum(relative_R.*relative_R); %�������
                term1 = R^(-6);
                term2 = R^(-3);
                
                if R < 0.25*L*L %С�ڽضϾ���
                    F = -24 * (2/R*term1 - 1/R*term2) .* relative_R; %������
                    %������ٶ�
                    for l = 1:3
                        storage_a(id_Itr,i,l) = storage_a(id_Itr,i,l) - F(l) ; %����i�ļ��ٶ�
                        storage_a(id_Itr,j,l) = storage_a(id_Itr,j,l) + F(l);  %����j�ļ��ٶ�
                    end
                    
                    %��������
                    storage_Potential(id_Itr,1) = storage_Potential(id_Itr,1) ...
                        + 4*(term1 - term2);
                    break
                else
                    continue
                end
            end
            
        end
    end
end

for i = 1:N
    
    %������ٶ�
    a = zeros(3,1);
    for l = 1:3
        a(l,1) = storage_a(id_Itr,i,l);
    end
    
    %����λ��
    sys(i,1).pos = sys(i,1).pos + time_step.*sys(i,1).v ...
        + 0.5*time_step*time_step.*a;
    
    %Խ����
    if sys(i,1).pos(1,1)>L || sys(i,1).pos(1,1)<0
        sys(i,1).pos(1,1) = mod(sys(i,1).pos(1,1),L);
    end
    if sys(i,1).pos(2,1)>L || sys(i,1).pos(2,1)<0
        sys(i,1).pos(2,1) = mod(sys(i,1).pos(2,1),L);
    end
    if sys(i,1).pos(3,1)>L || sys(i,1).pos(3,1)<0
        sys(i,1).pos(3,1) = mod(sys(i,1).pos(3,1),L);
    end
end
    
%���¼��ٶ�
if id_Itr ~= n_Itr
for i = 1:N
    for j = i+1:N
        for k = 1 : N  %����������
            relative_R = sys(i,1).pos - sys(j,1).pos + direc_vec(:,k);
            R = sum(relative_R.*relative_R); %�������
            term1 = R^(-6);
            term2 = R^(-3);
            
            if R < 0.25*L*L %С�ڽضϾ���
                F = -24 * (2/R*term1 - 1/R*term2) .* relative_R; %������
                %������ٶ�
                for l = 1:3
                    storage_a(id_Itr+1,i,l) = storage_a(id_Itr+1,i,l) - F(l) ; %����i�ļ��ٶ�
                    storage_a(id_Itr+1,j,l) = storage_a(id_Itr+1,j,l) + F(l);  %����j�ļ��ٶ�
                end
                
                %��������
                storage_Potential(id_Itr+1,1) = storage_Potential(id_Itr+1,1) ...
                    + 4*(term1 - term2);
                break
            else
                continue
            end
        end
        
    end
end
end

%������ٶ�   
for i = 1:N
    a = zeros(3,1);
    if id_Itr==n_Itr
        for l = 1:3
            a(l,1) = storage_a(id_Itr,i,l);
        end
    else
        for l = 1:3
            a(l,1) = storage_a(id_Itr,i,l) + storage_a(id_Itr+1,i,l);
        end
    end
    %�����ٶ� 
    sys(i,1).v = sys(i,1).v + 0.5*time_step.*a;    
    
    %��������
    storage_Ek(id_Itr,1) = storage_Ek(id_Itr,1) + sys(i,1).Ek;
    
    
end
%�����¶�
storage_temp(id_Itr,1) = 2*storage_Ek(id_Itr,1)/3/N;
end