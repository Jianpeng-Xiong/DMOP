classdef JY6 < PROBLEM
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
        function obj = JY6()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
%% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
           A = 0.1;
           W = 3;
           G=sin(0.5*pi*obj.Global.t);
           K = 2*floor(10*abs(G));
           y=PopDec(:,2:end)-G;
           g = 1+sum((4*y(:,2:end).^2-cos(K*pi*y(:,2:end))+1),2);
           PopObj(:,1) =g.*( PopDec(:,1) + A*sin(W*pi.*PopDec(:,1)));
           PopObj(:,2) = g.*(1-PopDec(:,1)+A*sin(W*pi.*PopDec(:,1)));
        end
        %% Sample reference points on Pareto front
        %此处将alfa和belta设置为1，可改为0-2任意值，但二者必须相同
        function P = PF(obj,N)
                N=10000;
                x=(0:1/(N-1):1)';
                f(:,1)=x+0.1*sin(3*pi*x);
                f(:,2)=1-x+0.1*sin(3*pi*x);
                P = f;
        end
    end
end
