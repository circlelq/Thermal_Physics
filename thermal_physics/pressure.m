function pressure(id_Itr)
global N L sys direc_vec storage_temp  storage_pressure

%PRESSURE Summary of this function goes here
new_pressure = 0;

for i = 1:N-1
	for j = i+1:N
	    for k = 1:27
	        relative_R = sys(i,1).pos - sys(j,1).pos + direc_vec(:,k);
             
	        if ((relative_R(1,1)<0.1425*L*L) && (relative_R(2,1)<0.1425*L*L)...
                    && (relative_R(2,1)<0.1425*L*L))
                %С�ڽضϾ��룬0.1425=0.57 * (0.5 * 0.5)
                R = sum(relative_R.*relative_R); %计算距离
                term1 = R^(-6);
                term2 = R^(-3);
	        	F = -24 * (2/R*term1 - 1/R*term2) .* relative_R; %计算�?

	            %计算压强
	            new_pressure = new_pressure + dot(relative_R,F);
	            break
	        else
	            continue
	        end
	    end
	end
end

storage_pressure(id_Itr) = new_pressure;

end

