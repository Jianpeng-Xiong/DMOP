classdef UDF3< PROBLEM
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
        function obj = UDF3()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D =10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = [0 -1*ones(1,obj.Global.D-1)];
            obj.Global.upper    = [1 ones(1,obj.Global.D-1)];
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
          function PopObj = CalObj(obj, PopDec)
            G=sin(0.5*pi*obj.Global.t);
            N=size(PopDec,1);
            t=(1/2*obj.Global.D+0.1)*(sin(2*10*pi*PopDec(:,1))-2*obj.Global.D*abs(G));
            tmp=max([zeros(N,1) t],[],2);
            odd=0;even=0;a1=0;b1=1;a2=0;b2=1;
            for i = 2:obj.Global.D
                if mod(i,2)
                   y1=PopDec(:,i)-sin(6*pi*PopDec(:,1)+(i*pi)/obj.Global.D);
                   a1=a1+2*(y1).^2;
                   b1=b1.*cos((20*pi*y1)/sqrt(i));
                   odd=odd+1;
                else
                   y2=PopDec(:,i)-sin(6*pi*PopDec(:,1)+(i*pi)/obj.Global.D);
                   a2=a2+2*(y2).^2;
                   b2=b2.*cos(20*pi*y2/sqrt(i));
                   even=even+1;  
                end
            end
            PopObj(:,1) = PopDec(:,1)+tmp+(2*(4*a1-2*b1+2).^2)/odd;
            PopObj(:,2) =1-PopDec(:,1)+tmp+(2*(4*a2-2*b2+2).^2)/even;
        end

%         function PopObj = CalObj(obj, PopDec)
%             G=sin(0.5*pi*obj.Global.t);
%             N=size(PopDec,1);
%             t=(1/2*obj.Global.D+0.1)*(sin(2*obj.Global.D*pi*PopDec(:,1))-2*obj.Global.D*abs(G));
%             tmp=max([zeros(N,1) t],[],2);
%             odd=0;even=0;a1=0;b1=1;a2=0;b2=1;
%             for i = 2:obj.Global.D
%                 if mod(i,2)
%                    y1=PopDec(:,i)-sin(6*pi*PopDec(:,1)+(i*pi)/obj.Global.D);
%                    a1=a1+2*(y1).^2;
%                    b1=b1.*cos((20*pi*y1)/sqrt(i));
%                    odd=odd+1;
%                 else
%                    y2=PopDec(:,i)-sin(6*pi*PopDec(:,1)+(i*pi)/obj.Global.D);
%                    a2=a2+2*(y2).^2;
%                    b2=b2.*cos(20*pi*y2/sqrt(i));
%                    even=even+1;  
%                 end
%             end
%             PopObj(:,1) = 1-PopDec(:,1)+tmp+(2*(4*a1-2*b1+2).^2)/even;
%             PopObj(:,2) =PopDec(:,1)+tmp+(2*(4*a2-2*b2+2).^2)/odd;
%    end




        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%           P = obj.load();
             N=10000;
            t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT);
            G=sin(0.5*pi*t);
            for i=1:N
               P(i,1) =(2*i-1)/(2*N)+abs(G);
               P(i,2) =(2*i)/(2*N)+abs(G);
            end
              T= P(:,2);
              P(:,2)=T(end:-1:1);
        end
    end
end