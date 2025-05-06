classdef DF2 < PROBLEM
% <problem> <DF>
% The 14 test functions are for cec2018 competition o dynamic multiobjective optimisation.
% tauT --- 10 --- the number of generation for which t remains fixed
% nT --- 10 --- the number of distinct steps in t

%------------------------------- Reference --------------------------------
% The 14 test functions are for cec2018 competition o dynamic multiobjective optimisation.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    properties(Access = public)
        tauT;   % the number of generation for which t remains fixed
        nT;  % the number of distinct steps in t
    end
    methods
        %% Initialization 
        function obj = DF2()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            G = abs(sin(0.5*pi*obj.Global.t));
            r=1+floor((obj.Global.D-1)*G);
            tmp=setdiff(1:obj.Global.D,r);
            a=length(tmp);g=0;
            for i=1:a
              g=g+(PopDec(:,tmp(i))-G).^2;
            end
            g=1+g;
            PopObj(:,1) = PopDec(:,r);
            PopObj(:,2) = g.*(1-(PopDec(:,r)./g).^0.5);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%             P = obj.load();
                N=10000;
               x=linspace(0,1,N);
               P(:,1)=x;
               P(:,2)=1-P(:,1).^0.5;
        end
    end
end