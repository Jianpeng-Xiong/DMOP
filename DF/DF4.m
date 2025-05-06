classdef DF4 < PROBLEM
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
        function obj = DF4()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = -2*ones(1,obj.Global.D);
            obj.Global.upper    = 2*ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            a = sin(0.5*pi*obj.Global.t);
            b=1+abs(cos(0.5*pi*obj.Global.t));
            c=max(abs(a),a+b);
            H=1.5+a;
            t=0;
            for i=2:obj.Global.D
                t=t+(PopDec(:,i)-(a*PopDec(:,1).^2)/(i*c^2)).^2;
            end
            g=1+t;
            PopObj(:,1) = (g.*abs(PopDec(:,1)-a).^H)+2*obj.Global.t; %注意：加了+2t
            PopObj(:,2) = g.*abs(PopDec(:,1)-a-b).^H+2*obj.Global.t; %注意：加了+2t
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
%             P = obj.load();
             N=10000;
             t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
             a=sin(0.5*pi*t);
             b=1+abs(cos(0.5*pi*t));     
             x=linspace(a,a+b,N);
             H=1.5+a;
             P(:,1)=abs(x-a).^H+2*t;
             P(:,2)=abs(x-a-b).^H+2*t;
        end
    end
end