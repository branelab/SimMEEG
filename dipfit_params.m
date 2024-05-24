classdef dipfit_params < handle
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        fitType             (1,1) string {mustBeMember(fitType, ["LeadField" "Continuous" "Multistage"])} = "LeadField";
        spatiotemp          (1,1) logical {mustBeMember(spatiotemp, [0 1])} = 1;   % (true)=spatiotemporal fitting across data samples (false)=moving dipole fitting across data samples;
        optimType           (1,1) string {mustBeMember(optimType, ["fminsearch", "fmincon", "fminmax", "Genetic", "GeneticMulti", "ParticleSwarm", "Pareto", "Pattern", "Surrogate", "SimulatedAnnealing"])} = "fminsearch";
        maxIter             (1,1) double {mustBeInRange(maxIter, 0, 1e10)} = 1000;  % max iterations for optimizations e.g., [1000];
        maxTime             (1,1) double {mustBeInRange(maxTime, 0, 1e10)} = 1000;  % max iterations for optimizations e.g., [1000];
        parallel_compute    (1,1) logical {mustBeMember(parallel_compute, [0 1])} = 1;  
        verbose             (1,1) string {mustBeMember(verbose, ["none",  "off",  "iter",  "diagnose",  "final"])} = "iter";  
        step_tolerance      (1,1) double {mustBeInRange(step_tolerance, 0, 1e10)} = .05;  % max iterations for optimizations e.g., [1000];
        fcn_tolerance       (1,1) double {mustBeInRange(fcn_tolerance, 0, 1e10)} = .05;  % max iterations for optimizations e.g., [1000];
    end

    methods
        %% Instantiate "dipole" with inputs from user - else it creates default 
        function obj = dipfit_params(spatiotemp, optimType, maxIter, verbose)
            if nargin==4    
                obj.spatiotemp = spatiotemp; 
                obj.optimType = optimType;
                obj.maxIter = maxIter;
                obj.verbose = verbose;
            end
        end

    end
end