classdef Particle
    properties
        pos = zeros(3,1);
        v = zeros(3,1);
        Ek = 0;
    end
    
    methods
        function self = set.pos(self,position)
        %position: 列向量
        %function: 设置点的pos属性
        self.pos(1,1) = position(1,1);
        self.pos(2,1) = position(2,1);
        self.pos(3,1) = position(3,1);
        end
         
        function self = set.Ek(self,vel)
        %vel: 列向量
        %function: 设置动能
        self.Ek = 0.5 * sum(vel.*vel);
        end
        
        function self = setEk(self,vel)
        self.Ek = vel;
        end
        
        function self = set.v(self,velocity)
        %velocity: 列向量
        %设置粒子的速度
        self.v(1,1) = velocity(1,1);
        self.v(2,1) = velocity(2,1);
        self.v(3,1) = velocity(3,1);
        self = self.setEk(velocity);  %更新速度后自动更新动能
        end
    end
end