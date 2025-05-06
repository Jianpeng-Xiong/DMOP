classdef JY10 < PROBLEM
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
        function obj = JY10()
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
        % rou represents the frequency of type change
        % A and W are two parameters to control the local shape of the POF
        % alfa and belta are parameters that control the overall shape of the POF
        function PopObj = CalObj(obj, PopDec)
            A=0.5;
            W=6;
            rou=5;
            temp=0;
            G=abs(sin(0.5*pi*obj.Global.t));
            sigma=mod((floor(obj.Global.gen/(obj.tauT*rou))+randperm(3,1)),3);
            alfa=1+sigma*G;
            belta=1+sigma*G;
            for i = 2:obj.Global.D
                 temp = temp + (PopDec(:,i)+sigma - G).^2;
            end
            g=temp;
            PopObj(:,1) =((1+g).*((PopDec(:,1)+A*sin(W*pi*PopDec(:,1))).^alfa))';
            PopObj(:,2) =((1+g).*((1-PopDec(:,1)+A*sin(W*pi*PopDec(:,1))).^belta))';
        end
    %% Sample reference points on Pareto front
        function P = PF(obj,~)
                 A=0.5;
                 W=6;
                 rou=5;
                 N=1000;
                 sigma=mod((floor(obj.Global.gen/(obj.tauT*rou))+randperm(3,1)),3);
                 G=abs(sin(0.5*pi*obj.Global.t));
                 alfa=1+sigma*G;
                 belta=1+sigma*G;
                 gt=(Global.D-1).*((1-G).^2);
                 P(:,1) = (0:1/(N-1):1)';
                 f1=((1+g).*(P(:,1)+A*sin(W*pi*P(:,1)))^alfa);
                 f2=((1+g).*(1-P(:,1)+A*sin(W*pi*P(:,1)))^belta);
                 P(:,2) =(1+gt).*(1+2*A*sin(W*pi*((f1^(1/alfa)-f2^(1/belta))/(2*(1+gt))+0.5)));
        end
    end
end