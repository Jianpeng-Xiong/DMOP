classdef FDA4 < PROBLEM
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
        function obj = FDA4()
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
          function f=CosX(x)
            temp=sin(x);
            f=sqrt(1-temp.*temp);
          end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            G=abs(sin(0.5*pi*obj.Global.t));
            temp = 0;
            for i = obj.Global.M:obj.Global.D
                temp = temp + (PopDec(:,i) - G).^2;
            end
            g =temp;
            m=obj.Global.M;
            for j=1:obj.Global.M
               PopObj(:,j)=1+g;
               for k=1:m-1
                   PopObj(:,j)=PopObj(:,j).*sqrt(1-sin((PopDec(:,k)*pi*0.5)).*sin((PopDec(:,k)*pi*0.5)));
               end
               if j>1
                  PopObj(:,j)= PopObj(:,j).*sin(0.5*pi*PopDec(:,m));
               end
               m=m-1;
            end
%             PopObj(:,1) = (1+g).*cos(0.5*pi*PopDec(:,2)).*cos(0.5*pi*PopDec(:,1));
%             PopObj(:,2) = (1+g).*cos(0.5*pi*PopDec(:,2)).*sin(0.5*pi*PopDec(:,1));
%             PopObj(:,3) = (1+g).*sin(0.5*pi*PopDec(:,2));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%             P = obj.load();  
               N=10050;
              [t,~]=UniformPoint(N,3);
              index=randperm(size(t,1),10000);
              t=t(index,:);
              P=t./repmat(sqrt(sum(t.^2,2)),1,3);
        end
    end
end