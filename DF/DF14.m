classdef DF14 < PROBLEM
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
        tauT;     % the number of generation for which t remains fixed  
        nT;       % the number of distinct steps in t  
    end
    methods
        %% Initialization 
        function obj = DF14()
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = [zeros(1,2) -1*ones(1,obj.Global.D-2)];
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            G=sin(0.5*pi*obj.Global.t);
            g=1+sum((PopDec(:,3:end)-G).^2,2);
            y=0.5+G*(PopDec(:,1)-0.5);
            PopObj(:,1) =g.*(1-y+0.05*sin(6*pi*y));
            PopObj(:,2) =g.*(1-PopDec(:,2)+0.05*sin(6*pi*PopDec(:,2))).*(y+0.05*sin(6*pi*y));
            PopObj(:,3) =g.*(PopDec(:,2)+0.05*sin(6*pi*PopDec(:,2))).*(y+0.05*sin(6*pi*y));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%             P = obj.load();
               g=1;N=15000;
                t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
               [x1,x2]=meshgrid(linspace(0,1,N/150));
               G=sin(0.5*pi*t);
               y=0.5+G*(x1-0.5);
               f1=g.*(1-y+0.05*sin(6*pi*y));
               f2=g.*(1-x2+0.05*sin(6*pi*x2)).*(y+0.05*sin(6*pi*y));
               f3=g.*(x2+0.05*sin(6*pi*x2)).*(y+0.05*sin(6*pi*y));
               [h]=get_PF({f1,f2, f3}, false);
               P=h;
        end
    end
end