classdef F10< PROBLEM
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
        function obj = F10()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D)*5;
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            a=2*cos(obj.Global.t*pi)+2;
            b=2*sin(2*pi*obj.Global.t)+2;
            H = 1.25+0.75*sin(pi*obj.Global.t);
            temp1=0;temp2=0;
            for i = 1:obj.Global.D
                if mod(i,2)
                   temp1 = temp1 + (PopDec(:,i)-b-(abs(PopDec(:,1)-a)).^(H+i/obj.Global.D)).^2;
                else 
                   temp2 = temp2 + (PopDec(:,i)-b-1+(abs(PopDec(:,1)-a)).^(H+i/obj.Global.D)).^2;
                end
            end
            PopObj(:,1) = abs(PopDec(:,1)-a).^H+temp1;
            PopObj(:,2) = abs(PopDec(:,1)-a-1).^H+temp2;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
              N=10000;
              t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
              H=1.25+0.75*sin(pi*t);
              s=linspace(0,1,N)';
              P(:,1)=s.^H;
              P(:,2) =(1-s).^H;
        end
    end
end