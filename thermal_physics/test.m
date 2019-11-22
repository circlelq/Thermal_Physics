for i = 1:27
    plot3(sys(i,1).v(1,1),sys(i,1).v(2,1),sys(i,1).v(3,1),'o');
    hold on
end