classdef JY8 < PROBLEM
% <problem> <JY>
% Benchmark DMOP proposed by Shouyong Jiang and Shengxiang Yang
% tauT --- 10 --- the number of generation for which t remains fixed
% nT --- 10 --- the number of distinct steps in t
%------------------------------- Reference --------------------------------
% Evolutionary Dynamic Multiobjective Optimization: Benchmarks and 
% Algorithm Comparisons
% Shouyong Jiang and Shengxiang Yang, Senior Member, IEEE，2015
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%------------------------------- Compile --------------------------------
% This dynamic multi-objective optimization benchmark function(JY) is
% written by Jianpeng Xiong in 2023.
% This code is for reference only and cannot be directly used on 
% publicly available platforms.
%--------------------------------------------------------------------------
    properties(Access = public)
        tauT;  % the number of generation for which t remains fixed
        nT; % the number of distinct steps in t
    end
    methods
        %% Initialization
        function obj = JY8()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = [0 -1*ones(1,obj.Global.D-1)];
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            A=0.05;
            W=6;
            G = sin(0.5 * pi * obj.Global.t);
            b=10-9.8*abs(G);
            a=2/b;
            g = 1+sum(PopDec(:,2:end).^2,2);
            PopObj(:,1) =g.*( PopDec(:,1) + A*sin(W*pi.*PopDec(:,1))).^a;
            PopObj(:,2) = g.*(1-PopDec(:,1)+A*sin(W*pi.*PopDec(:,1))).^b;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            N=10000;
            G = sin(0.5 * pi * obj.Global.t);
            b=10-9.8*abs(G);
            a=2/b;
            x=(0:1/(N-1):1)';
            f(:,1)=x+0.05*sin(6*pi*x);
            f(:,2)=1-x+0.05*sin(6*pi*x);
            P = f;
        end
    end
end
                 
            
