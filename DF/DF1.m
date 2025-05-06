classdef DF1 < PROBLEM
% <problem> <DF>
% The 14 test functions are for cec2018 competition on dynamic multiobjective optimisation.
% tauT --- 10 --- the number of generation for which t remains fixed
% nT --- 10 --- the number of distinct steps in t

%------------------------------- Reference --------------------------------
% Jiang S, Yang S. Evolutionary Dynamic Multiobjective Optimization: 
% Benchmarks and Algorithm Comparisons. IEEE Trans Cybern. 2017 Jan;
% 47(1):198-211. doi: 10.1109/TCYB.2015.2510698. Epub 2016 Jan 13. PMID: 26766387.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
% This dynamic multi-objective optimization benchmark function(DF) is
% written by Jianpeng Xiong in 2023.
% This code is for reference only and cannot be directly used on 
% publicly available platforms.
%--------------------------------------------------------------------------
    properties(Access = public)
        tauT;   % the number of generation for which t remains fixed
        nT;  % the number of distinct steps in t
    end
    methods
        %% Initialization 
        function obj = DF1()
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
        function PopObj = CalObj(obj,PopDec)
            G = abs(sin(0.5*pi*obj.Global.t));
            H = 0.75*sin(0.5*pi*obj.Global.t)+1.25;
            g = 1 + sum((PopDec(:,2:end)-G).^2, 2);
            PopObj(:,1) = PopDec(:,1);
            PopObj(:,2) = g.*(1-(PopDec(:,1)./g).^H);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj,N)
            N = 10000;
            t= 1/obj.nT * fix((obj.Global.gen-50-obj.tauT)/obj.tauT); 
            H = 0.75*sin(0.5*pi*t)+1.25;
            P(:,1) = (0:1/(N-1):1)';
            P(:,2) = 1 - P(:,1).^H;
        end
    end
end