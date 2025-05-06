classdef FDA2 < PROBLEM
% <problem> <FDA>
% Benchmark MOP proposed by Farina et al.
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
        function obj = FDA2()
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
      PopObj(:,1) = PopDec(:,1);
      g = 1+sum(PopDec(:,2:obj.Global.D-8).^2,2);
      H = 2*sin(0.5*pi*(obj.Global.t-1));
      h = 1-(PopObj(:,1)./g).^(2.^(H+sum((PopDec(:,obj.Global.D-7:end)-H/4).^2,2)));
      PopObj(:,2)=g.*h;     
     end
     
    %% Sample reference points on Pareto front
        function P = PF(obj,N)
            N=10000;
            % 计算PS解集
            PS = zeros(N,obj.Global.D);
            t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
            PS(:,1) = linspace(0,1,N);
            H = 2*sin(0.5*pi*(t-1));
            PS(:,obj.Global.D-7:end) = repmat(H/4,N,8);  
            % 计算PS计算PF
            P(:,1) = PS(:,1);
            g = 1+sum(PS(:,2:obj.Global.D-8).^2,2);
            h = 1-(P(:,1)./g).^(2.^(H+sum((PS(:,obj.Global.D-7:end)-H/4).^2,2)));
            P(:,2) =g.*h;
            P= rm_dominated(P);
        end
   
    end
end