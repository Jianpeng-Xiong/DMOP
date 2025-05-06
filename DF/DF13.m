classdef DF13 < PROBLEM
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
        function obj = DF13()
            obj.Global.M = 3;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = [0,0,-1*ones(1,obj.Global.D-2)];
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            Gt = sin(0.5*pi*obj.Global.t);
            p = floor(6*Gt);
            g = 1+sum((PopDec(:,3:end)-Gt).^2, 2);
            PopObj(:,1) = g.*cos(0.5*pi*PopDec(:,1)).^2;
            PopObj(:,2) = g.*cos(0.5*pi*PopDec(:,2)).^2;
            PopObj(:,3) = g.*(sin(0.5*pi*PopDec(:,1)).^2+sin(0.5*pi*PopDec(:,1)).*cos(p*pi*PopDec(:,1)).^2 +...
                sin(0.5*pi*PopDec(:,2)).^2+sin(0.5*pi*PopDec(:,2)).*cos(p*pi*PopDec(:,2)).^2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,~)
            t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT);
            %N = floor(sqrt(N));
            N = 15000;
            [x1,x2]=meshgrid(linspace(0,1,N/150));
            Gt = sin(0.5*pi*t);
            p = floor(6*Gt);
            [s1, s2] = size(x1);
            g = 1;
            f1 = reshape(g.*cos(0.5*pi*x1).^2, 1,s1*s2);
            f2 = reshape(g.*cos(0.5*pi*x2).^2, 1, s1*s2);
            f3 = reshape(g.*(sin(0.5*pi*x1).^2+sin(0.5*pi*x1).*cos(p*pi*x1).^2 +...
                sin(0.5*pi*x2).^2+sin(0.5*pi*x2).*cos(p*pi*x2).^2), 1, s1*s2); 
            [h]=get_PF({f1,f2, f3}, false);
            P=h;
        end
    end
end
