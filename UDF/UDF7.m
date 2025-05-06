classdef UDF7< PROBLEM
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
        function obj = UDF7()
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = [0 0 -2*ones(1,obj.Global.D-2)];
            obj.Global.upper    = [1 1 2*ones(1,obj.Global.D-2)];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            G=sin(0.5*pi*obj.Global.t);
            rt=1+abs(G);
            odd=0;trip=0;even=0;temp1=0;temp2=0;temp3=0;
            for i = 3:obj.Global.D
                y1=(PopDec(:,i)-2*PopDec(:,2).*sin(2*pi*PopDec(:,1)+i*pi/obj.Global.D)).^2;
                if mod(i,3)==0
                   temp1=temp1+y1;
                   odd=odd+1;
                elseif mod(i+1,3)==0
                   temp2=temp2+y1;
                   trip=trip+1;
                elseif mod(i-1,3)==0
                   temp3=temp3+y1;
                   even=even+1;
                end
            end
            PopObj(:,1) =rt*cos(0.5*pi*PopDec(:,1)).*cos(0.5*pi*PopDec(:,2))+2*temp1/odd+abs(G);
            PopObj(:,2) =rt*sin(0.5*pi*PopDec(:,1))+2*temp2/trip+abs(G);
            PopObj(:,3)=rt*cos(0.5*pi*PopDec(:,1)).*sin(0.5*pi*PopDec(:,2))+2*temp2/even+abs(G);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
             N=10050;
             t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT);
             G=sin(0.5*pi*t);
             R=1+abs(G);
             [F,~]=UniformPoint(N,3);
             index=randperm(size(F,1),10000);
             F=F(index,:);
             F=F./repmat(sqrt(sum(F.^2,2)),1,3);
             P=F.*R+G;
        end
    end
end