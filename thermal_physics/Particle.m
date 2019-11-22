classdef Particle
    properties
        pos = zeros(3,1);
        v = zeros(3,1);
        Ek = 0;
    end
    
    methods
        function self = set.pos(self,position)
        %position: ������
        %function: ���õ��pos����
        self.pos(1,1) = position(1,1);
        self.pos(2,1) = position(2,1);
        self.pos(3,1) = position(3,1);
        end
         
        function self = set.Ek(self,vel)
        %vel: ������
        %function: ���ö���
        self.Ek = 0.5 * sum(vel.*vel);
        end
        
        function self = setEk(self,vel)
        self.Ek = vel;
        end
        
        function self = set.v(self,velocity)
        %velocity: ������
        %�������ӵ��ٶ�
        self.v(1,1) = velocity(1,1);
        self.v(2,1) = velocity(2,1);
        self.v(3,1) = velocity(3,1);
        self = self.setEk(velocity);  %�����ٶȺ��Զ����¶���
        end
    end
end