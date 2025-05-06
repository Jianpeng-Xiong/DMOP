classdef F7< PROBLEM
% <problem> <F>
% Benchmark MOP proposed by Farina et al.
% tauT --- 5 --- the number of generation for which t remains fixed
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
    properties(Access = public)
        tauT;  % the number of generation for which t remains fixed
        nT; % the number of distinct steps in t
    end
    methods
        %% Initialization
        function obj = F7()
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
            a=1.7*(1-sin(pi*obj.Global.t))*sin(pi*obj.Global.t)+3.4;
            b=1.4*(1-sin(pi*obj.Global.t))*cos(pi*obj.Global.t)+2.1;
            H = 1.25+0.75*sin(pi*obj.Global.t);
            temp1=0;temp2=0;
            for i = 1:obj.Global.D
                if mod(i,2)
                   temp1 = temp1 + (PopDec(:,i)-b-1+(abs(PopDec(:,1)-a)).^(H+i/obj.Global.D)).^2;
                else 
                   temp2 = temp2 + (PopDec(:,i)-b-1+(abs(PopDec(:,1)-a)).^(H+i/obj.Global.D)).^2;
                end
            end
            PopObj(:,1) = abs(PopDec(:,1)-a).^H+temp1;
            PopObj(:,2) = abs(PopDec(:,1)-a-1).^H+temp2;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%              P = obj.load();
                N=10000;
                t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
                H=1.25+0.75*sin(pi*t);
                s=linspace(0,1,N)';
                P(:,1)=s.^H;
                P(:,2) =(1-s).^H;
        end
    end
end