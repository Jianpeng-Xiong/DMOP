classdef DF10 < PROBLEM
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
        function obj = DF10()
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
            H=2.25+2*cos(0.5*pi*obj.Global.t);
            tmp=sin(4*pi*(PopDec(:,1)+PopDec(:,2)))/(1+abs(G));
            g=1+sum((PopDec(:,3:end)-repmat(tmp,1,obj.Global.D-2)).^2,2);
            PopObj(:,1) =g.*sin(0.5*pi*PopDec(:,1)).^H;
            PopObj(:,2) =g.*sin(0.5*pi*PopDec(:,2)).^H.*cos(0.5*pi*PopDec(:,1)).^H;
            PopObj(:,3) =g.*cos(0.5*pi*PopDec(:,2)).^H.*cos(0.5*pi*PopDec(:,1)).^H;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
                    g=1;N=10000;
                   t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
                   [x1,x2]=meshgrid(linspace(0,1,N/100));
                   H=2.25+2*cos(0.5*pi*t);
                   f1=g*sin(0.5*pi*x1).^H;
                   f2=g*sin(0.5*pi*x2).^H.*cos(0.5*pi*x1).^H;
                   f3=g*cos(0.5*pi*x2).^H.*cos(0.5*pi*x1).^H;
                   [h]=get_PF({f1,f2, f3}, false);
                   P=h;                  
        end
    end
end