function plot3d()
%PLOT3D Summary of this function goes here
%   Detailed explanation goes here
global N L sys
for i = 1:N
	plot3(sys(i,1).pos(1,1),sys(i,1).pos(2,1),sys(i,1).pos(3,1),'o')
	hold on
end

