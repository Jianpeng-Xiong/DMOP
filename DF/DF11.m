classdef DF11 < PROBLEM
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
        function obj = DF11()
            obj.Global.M = 3;
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
            G=abs(sin(0.5*pi*obj.Global.t));
            g=1+G+sum((PopDec(:,3:end)-0.5*G*PopDec(:,1)).^2,2);
            y=pi*G/6+(pi/2-pi*G/3)*PopDec(:,1:2);
            PopObj(:,1) =g.*sin(y(:,1));
            PopObj(:,2) =g.*(sin(y(:,2)).*cos(y(:,1)));
            PopObj(:,3) =g.*(cos(y(:,2)).*cos(y(:,1)));
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
                N=100;
                t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
                G=abs(sin(0.5*pi*t));
                [X,Y] = meshgrid(0:1/(N-1):1, 0:1/(N-1):1);
                PS = [];
                for i = 1:size(X,1)
                    PS = [PS ;[X(:,i) Y(:,i)]];
                end
                for i = 1:size(PS,1)
                    for j = 3:obj.Global.D
                        PS(i,j) = 0.5*G*PS(i,1) ;
                    end
                end
               g=1+G+sum((PS(:,3:end)-repmat(0.5*G*PS(:,1),1,obj.Global.D-2)).^2,2);      
               y1=pi*G/6.0+(pi/2-pi*G/3.0)*PS(:,1);
               y2=pi*G/6.0+(pi/2-pi*G/3.0)*PS(:,2);
               f0=g.*sin(y1) ;
               f1=g.*sin(y2).*cos(y1);
               f2=g.*cos(y2).*cos(y1);
               P=[f0 f1 f2];
        end
    end
end