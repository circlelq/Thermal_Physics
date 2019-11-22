function Widom(id_Itr)
global N L sys direc_vec storage_temp chem_p boltz count

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
boltz = boltz + exp(-new_potential/storage_temp(id_Itr,1));
chem_p(count,1) = -storage_temp(id_Itr,1)*log(boltz/count);
end