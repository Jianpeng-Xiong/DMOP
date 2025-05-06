classdef DMOP2 < PROBLEM
% <problem> <DMOP>
% Benchmark MOP proposed by Chi-Keong Goh et al.
% tauT --- 10 --- the number of generation for which t remains fixed
% nT --- 10 --- the number of distinct steps in t

%------------------------------- Reference --------------------------------
% Chi-Keong Goh£¬ Kay Chen Tan. A Competitive-Cooperative Coevolutionary Paradigm
%for Dynamic Multiobjective Optimization.IEEE Transactions on Evolutionary Computation ( Volume: 13, Issue: 1, Feb. 2009)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%The code cannot run directly on your computer and is for reference only.

    properties(Access = public)
        tauT;  % the number of generation for which t remains fixed
        nT; % the number of distinct steps in t
    end
    methods
        %% Initialization
        function obj = DMOP2()
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
            G=sin(0.5*pi*obj.Global.t);
            H = 1.25+0.75*sin(0.5*pi*obj.Global.t);
            temp=0;
            for i = 2:obj.Global.D
                temp = temp + (PopDec(:,i)-G).^2;
            end
            g = 1 + temp;
            PopObj(:,1) = PopDec(:,1);
            h=1-(PopObj(:,1)./g).^H;
            PopObj(:,2) = g.*h;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%             P = obj.load();
             N=10000;
%              t= 1/obj.nT * fix((obj.Global.gen-obj.tauT)/obj.tauT);
            t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT);
            H=0.75*sin(0.5*pi*t)+1.25;
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1).^H; 
        end
    end
end