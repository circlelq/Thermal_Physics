function pressure(id_Itr)
global N L sys direc_vec storage_temp  storage_pressure

%PRESSURE Summary of this function goes here
new_pressure = 0;

for i = 1:N-1
	for j = i+1:N
	    for k = 1:27
	        relative_R = sys(i,1).pos - sys(j,1).pos + direc_vec(:,k);
	        R = sum(relative_R.*relative_R); %计算距离
	        term1 = R^(-6);
	        term2 = R^(-3);
	                
	        if R < 0.25*L*L %小于截断距离
	        	F = -24 * (2/R*term1 - 1/R*term2) .* relative_R; %计算力

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

