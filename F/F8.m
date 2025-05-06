classdef F8< PROBLEM
% <problem> <F>
% Benchmark MOP proposed by Aimin Zhou et al.
% tauT --- 10 --- the number of generation for which t remains fixed
% nT --- 10 --- the number of distinct steps in t

%------------------------------- Reference --------------------------------
% Aimin Zhou,Yaochu Jin,Qingfu Zhang. A Population Prediction Strategy for Evolutionary
% Dynamic Multiobjective Optimization.IEEE Transactions on Cybernetics ( Volume: 44, Issue: 1, Jan. 2014)
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
        function obj = F8()
             if isempty(obj.Global.M)
                obj.Global.M = 3;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = [zeros(1,2) -1*ones(1,obj.Global.D-2)];
            obj.Global.upper    = [ones(1,2) 2*ones(1,obj.Global.D-2)];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            H = 1.25+0.75*sin(pi*obj.Global.t);
            G=sin(0.5*pi*obj.Global.t);
            g=0;
            for i = 3:obj.Global.D
                g = g + (PopDec(:,i)-(((PopDec(:,1)+PopDec(:,2))/2).^H)-G).^2;
            end
            m=obj.Global.M;
            for j=1:obj.Global.M
                PopObj(:,j)=1+g;
                for k=1:m-1
                  PopObj(:,j)=  PopObj(:,j).*sqrt(1-sin(PopDec(:,k)*pi*0.5).*sin(PopDec(:,k)*pi*0.5));
                end
                if j>1
                   PopObj(:,j)=  PopObj(:,j).*sin(PopDec(:,m)*pi*0.5); 
                end
                m=m-1;
            end    

        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            N=10050;
            [t,~]=UniformPoint(N,3);
            index=randperm(size(t,1),10000);
            t=t(index,:);
            P=t./repmat(sqrt(sum(t.^2,2)),1,3);
        end
    end
end