classdef FDA3 < PROBLEM
% <problem> <FDA>
% Benchmark MOP proposed by Zitzler, Deb, and Thiele
% tauT --- 10 --- the number of generation for which t remains fixed
% nT --- 10 --- the number of distinct steps in t

%------------------------------- Reference --------------------------------
% Farina M, Deb K, Amato P. Dynamic multiobjective optimization problems: 
% test cases, approximations, and applications[J]. IEEE Transactions on 
% evolutionary computation, 2004, 8(5): 425-442.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
  properties(Access = public)
        tauT;  % the number of generation for which t remains fixed
         nT; % the number of distinct steps in t
    end
    methods
        %% Initialization
        function obj = FDA3()
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
        function PopObj = CalObj(obj,PopDec)
            F=10^(2*sin(0.5*pi*obj.Global.t));
            PopObj(:,1) = PopDec(:,1).^F;
            G=abs(sin(0.5*pi*obj.Global.t));
            temp = 0;
            for i = 2:obj.Global.D
                temp = temp + (PopDec(:,i) - G).^2;
            end
            g = 1 + G+temp;
            h = 1 - sqrt((PopObj(:,1)./g));
            PopObj(:,2) = g.*h;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%              P = obj.load(); 
               N=10000;
               t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT);
               G=abs(sin(0.5*pi*t));
               P(:,1)=(0:1/(N-1):1)';
               P(:,2)=(1+G)*(1-sqrt(P(:,1)/(1+G)));         
        end
    end
end