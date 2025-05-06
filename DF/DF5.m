classdef DF5 < PROBLEM
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
        function obj = DF5()
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
           w=floor(10*G);
            g=1+sum((PopDec(:,2:end)-G).^2,2);
            PopObj(:,1) = (g.*(PopDec(:,1)+0.02*sin(w*pi*PopDec(:,1))))+2*obj.Global.t; %注意：加了+2t
            PopObj(:,2) = (g.*(1-PopDec(:,1)+0.02*sin(w*pi*PopDec(:,1))))+2*obj.Global.t; %注意：加了+2t
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%             P = obj.load();
               N=10000;
              t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
              x=linspace(0,1,N);
              G=sin(0.5*pi*t);
              w=floor(10*G);
              P(:,1)=x+0.02*sin(w*pi*x)+2*t;
              P(:,2)=1-x+0.02*sin(w*pi*x)+2*t;
        end
    end
end