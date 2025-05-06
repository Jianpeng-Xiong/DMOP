classdef DF9 < PROBLEM
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
        function obj = DF9()
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            [obj.tauT, obj.nT] = obj.Global.ParameterSet(5, 10);
            obj.Global.lower    = [0 -1*ones(1,obj.Global.D-1)];
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj, PopDec)
            N=1+floor(10*abs(sin(0.5*pi*obj.Global.t)));
            g=1;
            for i=2:obj.Global.D
                tmp=PopDec(:,i)-cos(4*obj.Global.t+PopDec(:,1)+PopDec(:,i-1));
                g=g+tmp.^2;
            end
             PopObj(:,1) =g.*(PopDec(:,1)+max(0,(0.1+0.5/N)*sin(2*N*pi*PopDec(:,1))))+2*obj.Global.t;
             PopObj(:,2) =g.*(1-PopDec(:,1)+max(0,(0.1+0.5/N)*sin(2*N*pi*PopDec(:,1))))+2*obj.Global.t;
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
             N=20000;
             x=linspace(0,1,N);
             t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
             r=1+floor(10*abs(sin(0.5*pi*t)));
             f1=x+max(0, (0.1+0.5/r)*sin(2*r*pi*x))+2*t;
             f2=1-x+max(0, (0.1+0.5/r)*sin(2*r*pi*x))+2*t;
             [h]=get_PF({f1,f2}, false); 
             P(:,1)=h(2:end,1);
             P(:,2)=h(2:end,2);
        end   
    end
end