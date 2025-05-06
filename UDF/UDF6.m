classdef UDF6< PROBLEM
% <problem> <UDF>
% Benchmark MOP proposed by S Biswas et al.
% tauT --- 10 --- the number of generation for which t remains fixed
% nT --- 10 --- the number of distinct steps in t

%------------------------------- Reference --------------------------------
% S Biswas£¬S Das£¬PN Suganthan£¬CA Coello Coello. Evolutionary multiobjective optimization in dynamic environments: 
%A set of novel benchmark functions. 2014 IEEE Congress on Evolutionary Computation (CEC)
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% The code cannot run directly on your computer and is for reference only.

    properties(Access = public)
        tauT;  % the number of generation for which t remains fixed
        nT; % the number of distinct steps in t
    end
    methods
        %% Initialization
        function obj = UDF6()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = [0 -1*ones(1,obj.Global.D-1)];
            obj.Global.upper    = [1 ones(1,obj.Global.D-1)];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            G=sin(0.5*pi*obj.Global.t);
            mt=0.5+abs(G);
            tmp=(0.5/10+0.01)*abs(sin(2*10*pi*PopDec(:,1))-2*10*abs(G));
            odd=0;even=0;temp1=0;temp2=0;
            for i = 2:obj.Global.D
                if mod(i,2)
                   y1=PopDec(:,i)-sin(6*pi*PopDec(:,1)+(i*pi)/obj.Global.D);
                   temp1=temp1+(2*(y1).^2-cos(4*pi*y1)+1).^2;
                   odd=odd+1;
                else
                    y2=PopDec(:,i)-sin(6*pi*PopDec(:,1)+(i*pi)/obj.Global.D);
                    temp2=temp2+(2*(y2).^2-cos(4*pi*y2)+1).^2;
                   even=even+1;  
                end
            end
            PopObj(:,1) = PopDec(:,1)+tmp+(2*temp1)/odd;
            PopObj(:,2) =1-mt*PopDec(:,1)+tmp+(2*temp2)/even;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,~)
%              P = obj.load();
                N=10000;
               t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT);
               G=sin(0.5*pi*t);
               H=0.5+abs(G);
               for i=1:(2*N+1)
                   P(i,1) =(i-1)/(2*N)+abs(G);
                   P(i,2) =1-((i-1)/(2*N))*H+abs(G);
               end
        end
    end
end