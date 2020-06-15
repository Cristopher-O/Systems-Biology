%%  Genome-scale metabolic model (GEM) of Hydrogenophaga pseudoflava strain DSM 1084

%   The work herein presented was conducted in the systems biology group of the Department of Protein Science at the School of Engineering Sciences in Chemistry, Biotechnology,
%   and Health (CBH) KTH located at Science for Life Laboratory (SciLifeLab) in Solna, Stockholm (Sweden). It was the primary undertaking of a thesis project subsequently 
%   submitted for the degree of Master of Science in Engineering, Biotechnology at the Faculty of Engineering (LTH), Lund university in 2020. 

%   This thesis was titled "Reconstruction of Genome-Scale Metabolic Models with Concomitant Constraint-Based Modelling for Flux Prediction – a Case Study of Syngas Consuming 
%   Hydrogenophaga pseudoflava" and is available from Lund University Papers (LUP).

%   Author: Cristopher Ollagnier Widén
%   Supervisors: Paul Hudson (KTH), Ed van Niel (LTH), Kiyan Shabestary (KTH)

%%  Recommended software

%   Microsoft Windows 10 (64-bit)
%   MATLAB (v R2017b)                                                       
%   The MATLAB toolbox RAVEN (v 2.3.0) (Wang et al., 2018) 
%   The Gurobi Optimizer (v 8.0.1) solver. 
%   The libSBML MATLAB API (v 5.17.0) (Hucka et al., 2003)

%% Sequence-based GEM reconstruction

clc;                                                                        % Clears the Command Window
clear all;                                                                  % Clears the Workspace

%model=getKEGGModelForOrganism('hpse','','','',true,false,false);           % Allows for automatic sequence-based reconstruction of a GEM for H. pseudoflava utilizing the KEGG 
                                                                            % database. Spontaneous (not enzyme-catalyzed) reactions are included whereas reactions with
                                                                            % undefined stoichiometries are excluded. Reactions labelled as incomplete or general are
                                                                            % also excluded (type "help getKEGGModelForOrganism" for a full explanation).                                                                          
                                                                            % The resulting GEM is referred to as the first draft reconstruction (FDR).
                                                                            % 'hpse' is the KEGG Organism Code for H. psuedoflava strain DSM 1084, available at DDBJ/ENA/GenBank 
                                                                            % under the accessions CP037867 (chromosome) and CP037868 (megaplasmid pDSM1084). See also:
                                                                            % https://www.genome.jp/kegg-bin/show_organism?category=Hydrogenophaga%20pseudoflava&category_type=species
                                                     
load('HPseFDR.mat');                                                        % Loads a model generated from KEGG on the 10th of April 2020 using the command
                                                                            % "model=getKEGGModelForOrganism('hpse','','','',true,false,false);" which consisits of 1367 
                                                                            % reactions, 1583 metabolites and 914 genes.                                                                            
                                                                           
model.description=('Genome-scale metabolic model (GEM) of Hydrogenophaga pseudoflava strain DSM 1084');         % Adds a model description.
model.annotation.defaultLB=(-1001);                                                                             % Default lower bound is set to -1001 mmol gCDW^-1 h^-1.                        
model.annotation.defaultUB=(1001);                                                                              % Default upper bound is set to 1001 mmol gCDW^-1 h^-1.
model.annotation.givenName=('Cristopher');                                                                      % Adds the first name of the GEM author.
model.annotation.familyName=('Ollagnier Widén');                                                                % Adds the family name of the GEM author.
model.annotation.email=('c_widens@hotmail.com');                                                                % Adds the e-mail address to the GEM author.
model.annotation.organization=('Faculty of Engineering (LTH), Lund University');                                % Adds an organization.
model.annotation.note=('Genome-scale metabolic model (GEM) of Hydrogenophaga pseudoflava strain DSM 1084');     % Adds a note to the GEM user.
                                                                            
                                                                            % NB! This bacterial GEM is to be compartmentalized into two distinct compartments: 
                                                                            % i)    Cytosol (c), and 
                                                                            % ii)   Extracellular (e).
                                                                            
                                                                            % Constructing the two compartments: 
                                                                            
model.comps=[{'c'};{'e'}];
model.compNames=[{'Cytosol'};{'Extracellular'}];

                                                                            % All of the reactions in the FDR are initially assigned a lower bound of -1000 and an upper 
                                                                            % bound of 1000. These reaction bounds are later constrained in accordance with data on reaction
                                                                            % directionalities available from e.g. the MetaCyc database, where reaction directionalities have 
                                                                            % been curated manually by trained scientists (see below).

model.rev(:,1)=[1];                                                         % All reactions in the FDR are to be interpreted as reversible for now. This is done to harmonize
                                                                            % the system prior to the introduction of constraints pertaining to reaction directionalities 
                                                                            % which is done later.
                                                                            
model=setParam(model,'lb',model.rxns,[-1000]);                              % Sets the lower bound of all reactions in the FDR to -1000 mmol gCDW^-1 h^-1.
model=setParam(model,'ub',model.rxns,[1000]);                               % Sets the upper bound of all reactions in the FDR to 1000 mmol gCDW^-1 h^-1.

                                                        
                                                                            % All of the reactions in the FDR are initially assigned a confidence score of 4 to denote that
                                                                            % they were added to the model using sequence data as evidence for their existence. A confidence score
                                                                            % of 4 corresponds to a reaction whose reaction directionality have been manually curated.
                                                                            % Manual curation will be carried out below, and reactions whose reaction directionality couldn't
                                                                            % be constrained will instead be assigned a lower level confidence score of 3.                                                                           
                                                                            
model.rxnConfidenceScores=(randn((numel(model.rxns)),1));model.rxnConfidenceScores(:,1)=[4];                    % Assigns a confidence score of 4 to all reactions in the FDR.

FDRrcns=numel(model.rxns);                                                  % The amount of reactions currently assigned a confidence score is saved as 'FDRrcns' 
                                                                            % for later use.
                                                                            
                                                                            
                                                                            % Every reaction is assigned a confidence score Ci?{0, 1, 2, 3, 4, 5, 6, 7} to indicate the 
                                                                            % likelihood of its actual presence in the reactome of H. pseudoflava. The higher the 
                                                                            % confidence score, the better the evidence motivating its incorporation. Confidence scores
                                                                            % are assigned accordingly:
                                                                            
                                                                            % CRITERIA                                                                  CONFIDENCE SCORE                                                                            
                                                                            % Biochemical data (direct evidence)                                        7
                                                                            %   e.g. enzyme assays
                                                                            % Genetic data                                                              6
                                                                            %   e.g. knock-out/-in or overexpression analysis
                                                                            % Physiological data (indirect evidence)                                    5
                                                                            %   e.g. secretion products or defined medium requirements, 
                                                                            %   transport-, and exchange reactions
                                                                            % Sequence data (genome annotation)                                         4
                                                                            %   (reaction directionality successfully curated)
                                                                            % Sequence data (genome annotation)                                         3
                                                                            %   (reaction directionality not curated)
                                                                            % Modelling data – required for functional model, hypothetical reaction     2
                                                                            %   (more likely)
                                                                            % Modelling data – required for functional model, hypothetical reaction     1
                                                                            %   (less likely)
                                                                            % No evidence                                                               0
                                                                            %   e.g. fake reactions, lumped reactions, etc.
                                                                            

                                                                            % simplifyModel is run to deduce the amount of dead-end metabolites and dead-end reactions
                                                                            % in the FDR:

[reducedModel, deletedReactions, deletedMetabolites]=simplifyModel(model,false,false,false,true);

size(deletedMetabolites);                                                   % 1149 out of 1583 metabolites are 'orphan metabolites' - i.e., metabolites that can either 
                                                                            % only be produced or only consumed, but not both. In other words, dead-end metabolites.
                                                                            
size(deletedReactions);                                                     % 799 out of 1367 reactions are dead-end reactions.                                                                         
                                                                            

                                                                            % An extensive literature search was carried out with the objective to identify experimentally
                                                                            % verified sole carbon sources that H. pseudoflava is confirmed to grow on. These metabolites
                                                                            % were then categorized into:
                                                                            
                                                                            % (1)  those that alredy exist in the KEGG-based automatic first draft reconstruction (FDR), and
                                                                            % (2)  those that do not already exist in the FDR.
                                                                            
                                                                            % Adding transport and consuming exchange reactions for the carbon source metabolites in (1) is 
                                                                            % relatively straight-forward as their presence is already accomodated for in the FDR. The 
                                                                            % addition of the metabolites in (2) requires concomitant addition of intracellular reactions 
                                                                            % not captured by the sequence-based FDR in order to account for their uptake. It should be 
                                                                            % noted that performing reconstruction so as to ensure that the GEM encapsulates all of the 
                                                                            % experimentally verified sole carbon sources and the reactions whose presence they implicitly 
                                                                            % require are important in order to retain biological accuracy.   
                                                                            
                                                                            % NB! It may be possible for H. pseudoflava to excrete some of these metabolites which would
                                                                            % justify the addition of producing exchange reactions in addition to consuming exchange 
                                                                            % reactions. However, to begin with, solely consuming exchange reactions are created as 
                                                                            % their presence is confirmed in the literature. Unlike the availability of information 
                                                                            % giving supporting evidence regarding what sole carbon sources H. pseudoflava is
                                                                            % known to grow on, there is a lack of information regarding what metabolites it has the 
                                                                            % capacity to excrete - hence why only consuming exchange reactions are created first.
                                                                                                                                                        
                                                                            % To begin with, transport and consuming exchange reactions of confirmed sole carbon sources in 
                                                                            % category 1 are added in alphabetical order:                                                                           
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Acetate'},true,false);                  % Adds a reversible transport reaction for Acetate.
model.rxns(end,1)={'TRP_c->e_Acetate'};                                                     % Assigns an appropriate reaction ID to the transport reaction.
lastMetAcetate=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMetAcetate);                              % Adds a consuming exchange reaction for Acetate.  
model.rxns(end,1)={'EXC_IN_Acetate'};                                                       % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Alanine'},true,false);                % Adds a reversible transport reaction for D-Alanine.
model.rxns(end,1)={'TRP_c->e_D-Alanine'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Alanine.
model.rxns(end,1)={'EXC_IN_D-Alanine'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.

                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Alanine'},true,false);                % Adds a reversible transport reaction for L-Alanine.
model.rxns(end,1)={'TRP_c->e_L-Alanine'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Alanine.
model.rxns(end,1)={'EXC_IN_L-Alanine'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.

                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'4-Aminobutanoate'},true,false);         % Adds a reversible transport reaction for 4-Aminobutanoate.
model.rxns(end,1)={'TRP_c->e_4-Aminobutanoate'};                                            % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for 4-Aminobutanoate.
model.rxns(end,1)={'EXC_IN_4-Aminobutanoate'};                                              % Assigns an appropriate reaction ID to the exchange reaction.

                                                                        
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'5-Aminopentanoate'},true,false);        % Adds a reversible transport reaction for 5-Aminopentanoate.
model.rxns(end,1)={'TRP_c->e_5-Aminopentanoate'};                                           % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for 5-Aminopentanoate.
model.rxns(end,1)={'EXC_IN_5-Aminopentanoate'};                                             % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Arabitol'},true,false);               % Adds a reversible transport reaction for D-Arabitol.
model.rxns(end,1)={'TRP_c->e_D-Arabitol'};                                                  % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Arabitol.
model.rxns(end,1)={'EXC_IN_D-Arabitol'};                                                    % Assigns an appropriate reaction ID to the exchange reaction.

                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Arginine'},true,false);               % Adds a reversible transport reaction for D-Arginine.
model.rxns(end,1)={'TRP_c->e_D-Arginine'};                                                  % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Arginine.
model.rxns(end,1)={'EXC_IN_D-Arginine'};                                                    % Assigns an appropriate reaction ID to the exchange reaction.
                                                                            

[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Asparagine'},true,false);             % Adds a reversible transport reaction for L-Asparagine.
model.rxns(end,1)={'TRP_c->e_L-Asparagine'};                                                % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Asparagine.
model.rxns(end,1)={'EXC_IN_L-Asparagine'};                                                  % Assigns an appropriate reaction ID to the exchange reaction.

                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Aspartate'},true,false);              % Adds a reversible transport reaction for D-Aspartate.
model.rxns(end,1)={'TRP_c->e_D-Aspartate'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Aspartate.
model.rxns(end,1)={'EXC_IN_D-Aspartate'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.
                                                                            
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Citrate'},true,false);                  % Adds a reversible transport reaction for Citrate.
model.rxns(end,1)={'TRP_c->e_Citrate'};                                                     % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Citrate.
model.rxns(end,1)={'EXC_IN_Citrate'};                                                       % Assigns an appropriate reaction ID to the exchange reaction.

                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Cellobiose'},true,false);               % Adds a reversible transport reaction for Cellobiose.
model.rxns(end,1)={'TRP_c->e_Cellobiose'};                                                  % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Cellobiose.
model.rxns(end,1)={'EXC_IN_Cellobiose'};                                                    % Assigns an appropriate reaction ID to the exchange reaction.
                                                                           

[model, addedRxns]=addTransport(model,{'c'},{'e'},{'CO'},true,false);                       % Adds a reversible transport reaction for CO.
model.rxns(end,1)={'TRP_c->e_CO'};                                                          % Assigns an appropriate reaction ID to the transport reaction.
lastMetCO=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMetCO);                                   % Adds a consuming exchange reaction for CO.
model.rxns(end,1)={'EXC_IN_CO'};                                                            % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'CO2'},true,false);                      % Adds a reversible transport reaction for CO2.
model.rxns(end,1)={'TRP_c->e_CO2'};                                                         % Assigns an appropriate reaction ID to the transport reaction.
lastMetCO2=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMetCO2);                                  % Adds a consuming exchange reaction for CO2.
model.rxns(end,1)={'EXC_IN_CO2'};                                                           % Assigns an appropriate reaction ID to the exchange reaction.
                                                                            

[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Ethanolamine'},true,false);             % Adds a reversible transport reaction for Ethanolamine.
model.rxns(end,1)={'TRP_c->e_Ethanolamine'};                                                % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Ethanolamine.
model.rxns(end,1)={'EXC_IN_Ethanolamine'};                                                  % Assigns an appropriate reaction ID to the exchange reaction.
                                                                            
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Ethanol'},true,false);                  % Adds a reversible transport reaction for Ethanol.
model.rxns(end,1)={'TRP_c->e_Ethanol'};                                                     % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Ethanol.
model.rxns(end,1)={'EXC_IN_Ethanol'};                                                       % Assigns an appropriate reaction ID to the exchange reaction.

                                                    
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Fructose'},true,false);               % Adds a reversible transport reaction for D-Fructose.
model.rxns(end,1)={'TRP_c->e_D-Fructose'};                                                  % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Fructose.
model.rxns(end,1)={'EXC_IN_D-Fructose'};                                                    % Assigns an appropriate reaction ID to the exchange reaction.

                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Fumarate'},true,false);                 % Adds a reversible transport reaction for Fumarate.
model.rxns(end,1)={'TRP_c->e_Fumarate'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Fumarate.
model.rxns(end,1)={'EXC_IN_Fumarate'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.
                                                                         
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Galactose'},true,false);              % Adds a reversible transport reaction for D-Galactose.
model.rxns(end,1)={'TRP_c->e_D-Galactose'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Galactose.
model.rxns(end,1)={'EXC_IN_D-Galactose'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.
                     
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Gluconic acid'},true,false);          % Adds a reversible transport reaction for D-Gluconic acid.
model.rxns(end,1)={'TRP_c->e_D-Gluconic acid'};                                             % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Gluconic acid.      
model.rxns(end,1)={'EXC_IN_D-Gluconic acid'};                                               % Assigns an appropriate reaction ID to the exchange reaction.
                                                                                                                                                        
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Glucose'},true,false);                % Adds a reversible transport reaction for D-Glucose.
model.rxns(end,1)={'TRP_c->e_D-Glucose'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Glucose.
model.rxns(end,1)={'EXC_IN_D-Glucose'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.
                                                                           
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Glutamate'},true,false);              % Adds a reversible transport reaction for D-Glutamate.
model.rxns(end,1)={'TRP_c->e_D-Glutamate'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Glutamate.
model.rxns(end,1)={'EXC_IN_D-Glutamate'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.
                                                                            
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Glutamate'},true,false);              % Adds a reversible transport reaction for L-Glutamate.
model.rxns(end,1)={'TRP_c->e_L-Glutamate'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Glutamate.
model.rxns(end,1)={'EXC_IN_L-Glutamate'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Glycerate'},true,false);              % Adds a reversible transport reaction for D-Glycerate.
model.rxns(end,1)={'TRP_c->e_D-Glycerate'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Glycerate.            
model.rxns(end,1)={'EXC_IN_D-Glycerate'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Glycerol'},true,false);                 % Adds a reversible transport reaction for Glycerol.
model.rxns(end,1)={'TRP_c->e_Glycerol'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
lastMetGlycerol=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMetGlycerol);                             % Adds a consuming exchange reaction for Glycerol.
model.rxns(end,1)={'EXC_IN_Glycerol'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Glycolate'},true,false);                % Adds a reversible transport reaction for Glycolate.
model.rxns(end,1)={'TRP_c->e_Glycolate'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Glycolate.
model.rxns(end,1)={'EXC_IN_Glycolate'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Histidine'},true,false);              % Adds a reversible transport reaction for L-Histidine.
model.rxns(end,1)={'TRP_c->e_L-Histidine'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Histidine.
model.rxns(end,1)={'EXC_IN_L-Histidine'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'4-Hydroxybenzoate'},true,false);        % Adds a reversible transport reaction for 4-Hydroxybenzoate.
model.rxns(end,1)={'TRP_c->e_4-Hydroxybenzoate'};                                           % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for 4-Hydroxybenzoate.
model.rxns(end,1)={'EXC_IN_4-Hydroxybenzoate'};                                             % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'(R)-3-Hydroxybutanoate'},true,false);   % Adds a reversible transport reaction for (R)-3-Hydroxybutanoate.
model.rxns(end,1)={'TRP_c->e_(R)-3-Hydroxybutanoate'};                                      % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for (R)-3-Hydroxybutanoate.
model.rxns(end,1)={'EXC_IN_(R)-3-Hydroxybutanoate'};                                        % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Isoleucine'},true,false);             % Adds a reversible transport reaction for L-Isoleucine.
model.rxns(end,1)={'TRP_c->e_L-Isoleucine'};                                                % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Isoleucine.
model.rxns(end,1)={'EXC_IN_L-Isoleucine'};                                                  % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'(R)-Lactate'},true,false);              % Adds a reversible transport reaction for (R)-Lactate.
model.rxns(end,1)={'TRP_c->e_(R)-Lactate'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for (R)-Lactate.
model.rxns(end,1)={'EXC_IN_(R)-Lactate'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'(S)-Lactate'},true,false);              % Adds a reversible transport reaction for (S)-Lactate.
model.rxns(end,1)={'TRP_c->e_(S)-Lactate'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for (S)-Lactate.
model.rxns(end,1)={'EXC_IN_(S)-Lactate'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Leucine'},true,false);                % Adds a reversible transport reaction for L-Leucine.
model.rxns(end,1)={'TRP_c->e_L-Leucine'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Leucine.
model.rxns(end,1)={'EXC_IN_L-Leucine'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.

                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Lysine'},true,false);                 % Adds a reversible transport reaction for L-Lysine.
model.rxns(end,1)={'TRP_c->e_L-Lysine'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Lysine.
model.rxns(end,1)={'EXC_IN_L-Lysine'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'(R)-Malate'},true,false);               % Adds a reversible transport reaction for (R)-Malate.
model.rxns(end,1)={'TRP_c->e_(R)-Malate'};                                                  % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for (R)-Malate.
model.rxns(end,1)={'EXC_IN_(R)-Malate'};                                                    % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'(S)-Malate'},true,false);               % Adds a reversible transport reaction for (S)-Malate.
model.rxns(end,1)={'TRP_c->e_(S)-Malate'};                                                  % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for (S)-Malate. 
model.rxns(end,1)={'EXC_IN_(S)-Malate'};                                                    % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Maltose'},true,false);                  % Adds a reversible transport reaction for Maltose.
model.rxns(end,1)={'TRP_c->e_Maltose'};                                                     % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Maltose.
model.rxns(end,1)={'EXC_IN_Maltose'};                                                       % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Mannitol'},true,false);                 % Adds a reversible transport reaction for Mannitol.
model.rxns(end,1)={'TRP_c->e_Mannitol'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Mannitol. 
model.rxns(end,1)={'EXC_IN_Mannitol'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Ornithine'},true,false);              % Adds a reversible transport reaction for L-Ornithine.
model.rxns(end,1)={'TRP_c->e_L-Ornithine'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Ornithine. 
model.rxns(end,1)={'EXC_IN_L-Ornithine'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Phenylalanine'},true,false);          % Adds a reversible transport reaction for L-Phenylalanine.
model.rxns(end,1)={'TRP_c->e_L-Phenylalanine'};                                             % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Phenylalanine. 
model.rxns(end,1)={'EXC_IN_L-Phenylalanine'};                                               % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Proline'},true,false);                % Adds a reversible transport reaction for L-Proline.
model.rxns(end,1)={'TRP_c->e_L-Proline'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Proline.              
model.rxns(end,1)={'EXC_IN_L-Proline'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.
                                                                            
                                                                                                                                                      
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Putrescine'},true,false);               % Adds a reversible transport reaction for Putrescine.
model.rxns(end,1)={'TRP_c->e_Putrescine'};                                                  % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Putrescine.
model.rxns(end,1)={'EXC_IN_Putrescine'};                                                    % Assigns an appropriate reaction ID to the exchange reaction.
 

[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Pyruvate'},true,false);                 % Adds a reversible transport reaction for Pyruvate.
model.rxns(end,1)={'TRP_c->e_Pyruvate'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Pyruvate.
model.rxns(end,1)={'EXC_IN_Pyruvate'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Serine'},true,false);                 % Adds a reversible transport reaction for L-Serine.
model.rxns(end,1)={'TRP_c->e_L-Serine'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Serine.
model.rxns(end,1)={'EXC_IN_L-Serine'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Spermine'},true,false);                 % Adds a reversible transport reaction for Spermine.
model.rxns(end,1)={'TRP_c->e_Spermine'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Spermine.
model.rxns(end,1)={'EXC_IN_Spermine'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.

                                                                        
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Succinate'},true,false);                % Adds a reversible transport reaction for Succinate.
model.rxns(end,1)={'TRP_c->e_Succinate'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Succinate.     
model.rxns(end,1)={'EXC_IN_Succinate'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.
                                                                            
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Sucrose'},true,false);                  % Adds a reversible transport reaction for Sucrose.
model.rxns(end,1)={'TRP_c->e_Sucrose'};                                                     % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Sucrose.
model.rxns(end,1)={'EXC_IN_Sucrose'};                                                       % Assigns an appropriate reaction ID to the exchange reaction.


[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Sucrose (G00370)'},true,false);         % Adds a reversible transport reaction for Sucrose (G00370).
model.rxns(end,1)={'TRP_c->e_Sucrose (G00370)'};                                            % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Sucrose (G00370).  
model.rxns(end,1)={'EXC_IN_Sucrose (G00370)'};                                              % Assigns an appropriate reaction ID to the exchange reaction.
                                                                                                   
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Threonine'},true,false);              % Adds a reversible transport reaction for L-Threonine.
model.rxns(end,1)={'TRP_c->e_L-Threonine'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Threonine.
model.rxns(end,1)={'EXC_IN_L-Threonine'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.
                                                                                                                                                      
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Tryptophan'},true,false);             % Adds a reversible transport reaction for L-Tryptophan.
model.rxns(end,1)={'TRP_c->e_L-Tryptophan'};                                                % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Tryptophan.
model.rxns(end,1)={'EXC_IN_L-Tryptophan'};                                                  % Assigns an appropriate reaction ID to the exchange reaction.
                                                                              
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Tyrosine'},true,false);               % Adds a reversible transport reaction for L-Tyrosine.
model.rxns(end,1)={'TRP_c->e_L-Tyrosine'};                                                  % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Tyrosine.
model.rxns(end,1)={'EXC_IN_L-Tyrosine'};                                                    % Assigns an appropriate reaction ID to the exchange reaction.
                                                                             
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Urea'},true,false);                     % Adds a reversible transport reaction for Urea.
model.rxns(end,1)={'TRP_c->e_Urea'};                                                        % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Urea.
model.rxns(end,1)={'EXC_IN_Urea'};                                                          % Assigns an appropriate reaction ID to the exchange reaction.
                                                                            
            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Valine'},true,false);                 % Adds a reversible transport reaction for L-Valine.
model.rxns(end,1)={'TRP_c->e_L-Valine'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Valine.
model.rxns(end,1)={'EXC_IN_L-Valine'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.
                                                                            
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Xylose'},true,false);                 % Adds a reversible transport reaction for D-Xylose.
model.rxns(end,1)={'TRP_c->e_D-Xylose'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Xylose.    
model.rxns(end,1)={'EXC_IN_D-Xylose'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.


                                                                            % All transport and consuming exchange reactions for metabolites in category (1) are assigned 
                                                                            % a confidence score of 5. The fact that physiological data reported in the literature has 
                                                                            % provided indirect evidence supporting the presence of these reactions in H. pseudoflava 
                                                                            % justifies this relatively high score.

model.rxnConfidenceScores((FDRrcns+1):numel(model.rxns),:)=[5];


                                                                            % Incorporation of sole carbon sources belogning to category 2 - i.e., carbon sources which 
                                                                            % necessitates concomitant adding of non-annotated reaction(s) are performed next. KEGG was advised
                                                                            % for reaction associations of all the compounds in category 2 and one or sometimes several 
                                                                            % reactions were then added to the model. Reactions were added from the general, full KEGG model
                                                                            % downloaded below. Added non-annotated reactions were assigned a confidence score of 2 if a reaction 
                                                                            % was added manually when and if there was only one reaction associated with the metabolite in 
                                                                            % question in KEGG and all other compounds participating in that reaction was already 
                                                                            % accounted for in the FDR. Similarly, a lower confidence score of 1 were assigned to added 
                                                                            % non-annotated reactions in cases where a reaction was added manually amongst several possible 
                                                                            % alternative reactions. For a full justification for why a particular reaction was added, the reader 
                                                                            % is referred to the thesis. 
                                                                            
%keggModel=getModelFromKEGG([],false,false,false,false);                    % Downloads the full KEGG model.                                                              
load('keggModelCris.mat');                                                  % Loads the full KEGG model (downloaded 20200312).
keggModel.comps=[{'c'}];                                                    % Equalizing compartmentation nomenclature between models.
keggModel.compNames=[{'Cytosol'}];                                          % Equalizing compartmentation nomenclature between models.

                                                                            % Addition of sole carbon source Acetamide:
                                                                            
                                                                            % There is only one reaction associated with Acetamide in KEGG: R00321.
                                                                            
                                                                            % Adding reaction R00321 (Acetamide + H2O <=> Acetate + Ammonia) is sufficient to connect
                                                                            % Acetamide to the FDR.

model=addRxnsGenesMets(model,keggModel,{'R00321'},false,[],2);              % Adds the reaction R00321 (Confidence score: 2) along with the metabolite Acetamide.
model.rev((find(ismember(model.rxns,'R00321'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R00321'))),:)=[2];                     % Reaffirms the confidence score of 2. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Acetamide'},true,false);                % Adds a reversible transport reaction for Acetamide.
model.rxns(end,1)={'TRP_c->e_Acetamide'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_Acetamide'))),:)=[5];         % Assigns a confidence score of 5. 
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Acetamide.
model.rxns(end,1)={'EXC_IN_Acetamide'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_Acetamide'))),:)=[5];           % Assigns a confidence score of 5. 

                                                                            % Addition of sole carbon source L-Arabinose:
                                                                            
model=addRxnsGenesMets(model,keggModel,{'R02526'},false,[],1);              % Adds the reaction R02526 (Confidence score: 1). 
model.rev((find(ismember(model.rxns,'R02526'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R02526'))),:)=[1];                     % Reaffirms the confidence score of 1. 
model=addRxnsGenesMets(model,keggModel,{'R01757'},false,[],1);              % Adds the reaction R01757 (Confidence score: 1) along with the metabolite L-Arabinose.
model.rev((find(ismember(model.rxns,'R01757'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R01757'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'L-Arabinose'},true,false);              % Adds a reversible transport reaction for L-Arabinose.
model.rxns(end,1)={'TRP_c->e_L-Arabinose'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_L-Arabinose'))),:)=[5];       % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for L-Arabinose.
model.rxns(end,1)={'EXC_IN_L-Arabinose'};                                                   % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_L-Arabinose'))),:)=[5];         % Assigns a confidence score of 5.           
            
                                                                            % Addition of sole carbon source Butanoic acid:
                                                                            
model=addRxnsGenesMets(model,keggModel,{'R01176'},false,[],1);              % Adds the reaction R01176 (Confidence score: 1) along with the metabolite Butanoic acid.
model.rev((find(ismember(model.rxns,'R01176'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R01176'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Butanoic acid'},true,false);            % Adds a reversible transport reaction for Butanoic acid.
model.rxns(end,1)={'TRP_c->e_Butanoic acid'};                                               % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_Butanoic acid'))),:)=[5];     % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Butanoic acid.
model.rxns(end,1)={'EXC_IN_Butanoic acid'};                                                 % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_Butanoic acid'))),:)=[5];       % Assigns a confidence score of 5.

                                                                            % Addition of sole carbon source 1-Butanol:
                                                                            
model=addRxnsGenesMets(model,keggModel,{'R03544'},false,[],1);              % Adds the reaction R03544 (Confidence score: 1) along with the metabolite 1-Butanol.
model.rev((find(ismember(model.rxns,'R03544'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R03544'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'1-Butanol'},true,false);                % Adds a reversible transport reaction for 1-Butanol.
model.rxns(end,1)={'TRP_c->e_1-Butanol'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_1-Butanol'))),:)=[5];         % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for 1-Butanol.
model.rxns(end,1)={'EXC_IN_1-Butanol'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_1-Butanol'))),:)=[5];           % Assigns a confidence score of 5.


                                                                            % Addition of sole carbon source D-Glucosamine:
                                                                            
model=addRxnsGenesMets(model,keggModel,{'R01961'},false,[],1);              % Adds the reaction R01961 (Confidence score: 1) along with the metabolite D-Glucosamine.
model.rev((find(ismember(model.rxns,'R01961'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R01961'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Glucosamine'},true,false);            % Adds a reversible transport reaction for D-Glucosamine.
model.rxns(end,1)={'TRP_c->e_D-Glucosamine'};                                               % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_D-Glucosamine'))),:)=[5];     % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Glucosamine.
model.rxns(end,1)={'EXC_IN_D-Glucosamine'};                                                 % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_D-Glucosamine'))),:)=[5];       % Assigns a confidence score of 5.

                                                                            % Addition of sole carbon source 3-Hydroxybenzoate:
                                                                            
model=addRxnsGenesMets(model,keggModel,{'R01628'},false,[],1);              % Adds the reaction R01628 (Confidence score: 1).
model.rev((find(ismember(model.rxns,'R01628'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R01628'))),:)=[1];                     % Reaffirms the confidence score of 1. 
model=addRxnsGenesMets(model,keggModel,{'R02589'},false,[],1);              % Adds the reaction R02589 (Confidence score: 1) along with the metabolite 3-Hydroxybenzoate.
model.rev((find(ismember(model.rxns,'R02589'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R02589'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'3-Hydroxybenzoate'},true,false);        % Adds a reversible transport reaction for 3-Hydroxybenzoate.
model.rxns(end,1)={'TRP_c->e_3-Hydroxybenzoate'};                                           % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_3-Hydroxybenzoate'))),:)=[5]; % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for 3-Hydroxybenzoate.
model.rxns(end,1)={'EXC_IN_3-Hydroxybenzoate'};                                             % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_3-Hydroxybenzoate'))),:)=[5];   % Assigns a confidence score of 5.

                                                                            % Addition of sole carbon source Lactose:

model=addRxnsGenesMets(model,keggModel,{'R01100'},false,[],1);              % Adds the reaction R01100 (Confidence score: 1) along with the metabolite Lactose.
model.rev((find(ismember(model.rxns,'R01100'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R01100'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Lactose'},true,false);                  % Adds a reversible transport reaction for Lactose.
model.rxns(end,1)={'TRP_c->e_Lactose'};                                                     % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_Lactose'))),:)=[5];           % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Lactose.
model.rxns(end,1)={'EXC_IN_Lactose'};                                                       % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_Lactose'))),:)=[5];             % Assigns a confidence score of 5.


                                                                            % Addition of sole carbon source D-Lyxose:
                                                                            
                                                                            % There is only one reaction associated with L-Lyxose in KEGG: R01898.
                                                                            
                                                                            % Adding reaction R01898 (D-Xylulose <=> D-Lyxose) is sufficient to connect
                                                                            % D-Lyxose to the FDR.                                                                          
                                                                            
model=addRxnsGenesMets(model,keggModel,{'R01898'},false,[],2);              % Adds the reaction R01898 (Confidence score: 2) along with the metabolite D-Lyxose.
model.rev((find(ismember(model.rxns,'R01898'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R01898'))),:)=[2];                     % Reaffirms the confidence score of 2. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Lyxose'},true,false);                 % Adds a reversible transport reaction for D-Lyxose.
model.rxns(end,1)={'TRP_c->e_D-Lyxose'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_D-Lyxose'))),:)=[5];          % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Lyxose.
model.rxns(end,1)={'EXC_IN_D-Lyxose'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_D-Lyxose'))),:)=[5];            % Assigns a confidence score of 5.

                                                                            % Addition of sole carbon source D-Mannose:
                                                                            
model=addRxnsGenesMets(model,keggModel,{'R00877'},false,[],1);              % Adds the reaction R00877 (Confidence score: 1) along with the metabolite D-Mannose.
model.rev((find(ismember(model.rxns,'R00877'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R00877'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Mannose'},true,false);                % Adds a reversible transport reaction for D-Mannose.
model.rxns(end,1)={'TRP_c->e_D-Mannose'};                                                   % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_D-Mannose'))),:)=[5];         % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Mannose. 
model.rxns(end,1)={'EXC_IN_D-Mannose'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_D-Mannose'))),:)=[5];           % Assigns a confidence score of 5.

                                                                            % Addition of sole carbon source Propane-1-ol:

model=addRxnsGenesMets(model,keggModel,{'R02377'},false,[],1);              % Adds the reaction R02377 (Confidence score: 1) along with the metabolite Propane-1-ol.
model.rev((find(ismember(model.rxns,'R02377'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R02377'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Propane-1-ol'},true,false);             % Adds a reversible transport reaction for Propane-1-ol.
model.rxns(end,1)={'TRP_c->e_Propane-1-ol'};                                                % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_Propane-1-ol'))),:)=[5];      % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Propane-1-ol.
model.rxns(end,1)={'EXC_IN_Propane-1-ol'};                                                  % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_Propane-1-ol'))),:)=[5];        % Assigns a confidence score of 5.

                                                                            % Addition of sole carbon source D-Sorbitol:

model=addRxnsGenesMets(model,keggModel,{'R00875'},false,[],1);              % Adds the reaction R00875 (Confidence score: 1) along with the metabolite D-Sorbitol.
model.rev((find(ismember(model.rxns,'R00875'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R00875'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'D-Sorbitol'},true,false);               % Adds a reversible transport reaction for D-Sorbitol.
model.rxns(end,1)={'TRP_c->e_D-Sorbitol'};                                                  % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_D-Sorbitol'))),:)=[5];        % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for D-Sorbitol.  
model.rxns(end,1)={'EXC_IN_D-Sorbitol'};                                                    % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_D-Sorbitol'))),:)=[5];          % Assigns a confidence score of 5.
                                                                            
                                                                            % Addition of sole carbon source alpha,alpha-Trehalose:

model=addRxnsGenesMets(model,keggModel,{'R00010'},false,[],1);              % Adds the reaction R00010 (Confidence score: 1) along with the metabolite alpha,alpha-Trehalose.
model.rev((find(ismember(model.rxns,'R00010'))),1)=[1];                     % The reaction is to be interpreted as reversible for now.
model.rxnConfidenceScores((find(ismember(model.rxns,'R00010'))),:)=[1];                     % Reaffirms the confidence score of 1. 
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'alpha,alpha-Trehalose'},true,false);    % Adds a reversible transport reaction for alpha,alpha-Trehalose.
model.rxns(end,1)={'TRP_c->e_alpha,alpha-Trehalose'};                                       % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_alpha,alpha-Trehalose'))),:)=[5];         % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for alpha,alpha-Trehalose.
model.rxns(end,1)={'EXC_IN_alpha,alpha-Trehalose'};                                         % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_alpha,alpha-Trehalose'))),:)=[5];           % Assigns a confidence score of 5.


                                                                            % Finally, a consuming exchange reaction is added for hydrogen (H2), in order to account for 
                                                                            % the hydrogen-oxidizing activity which H. pseudoflava possesses:
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Hydrogen'},true,false);                 % Adds a reversible transport reaction for Hydrogen (H2).
model.rxns(end,1)={'TRP_c->e_Hydrogen'};                                                    % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_Hydrogen'))),:)=[5];          % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Hydrogen (H2). 
model.rxns(end,1)={'EXC_IN_Hydrogen'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_Hydrogen'))),:)=[5];            % Assigns a confidence score of 5.                                                                           
                                                                            

                                                                            % In preparation for various forthcoming tests, all of the consuming exchange reactions added 
                                                                            % above for the various possible sole carbon (and energy) sources are closed.
                                                                            
model=setParam(model,'eq',getExchangeRxns(model),0);                        % Closes all consuming and producing exchange reactions.

                                                                            % Being an obligate aerobe, H. Pseudoflava needs oxygen (O2) to survive. Therefore, a
                                                                            % transport reaction is added along with consuming and producing exchange reactions 
                                                                            % for oxygen:

[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Oxygen'},true,false);                   % Adds a reversible transport reaction for Oxygen (O2).
model.rxns(end,1)={'TRP_c->e_Oxygen'};                                                      % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_Oxygen'))),:)=[5];            % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Oxygen (O2).
model.rxns(end,1)={'EXC_IN_Oxygen'};                                                        % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_Oxygen'))),:)=[5];              % Assigns a confidence score of 5.
[model, addedRxns]=addExchangeRxns(model,'out',lastMet);                                    % Adds a producing exchange reaction for Oxygen (O2).
model.rxns(end,1)={'EXC_OUT_Oxygen'};                                                       % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_OUT_Oxygen'))),:)=[5];             % Assigns a confidence score of 5.

                                                                            % Likewise, H2O is needed and so a transport reaction for water is added along with
                                                                            % a consuming and producing exchange reaction:
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'H2O'},true,false);                      % Adds a reversible transport reaction for H2O.
model.rxns(end,1)={'TRP_c->e_H2O'};                                                         % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_H2O'))),:)=[5];               % Assigns a confidence score of 5.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for H2O.
model.rxns(end,1)={'EXC_IN_H2O'};                                                           % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_H2O'))),:)=[5];                 % Assigns a confidence score of 5.
[model, addedRxns]=addExchangeRxns(model,'out',lastMet);                                    % Adds a producing exchange reaction for H2O.      
model.rxns(end,1)={'EXC_OUT_H2O'};                                                          % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_OUT_H2O'))),:)=[5];                % Assigns a confidence score of 5.

                                                                            % Similarly, a nitrogen source in the form of Ammonia is incorporated into the GEM
                                                                            % through addition of a transport reaction along with a consuming exchange reaction:
                                                                            
[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Ammonia'},true,false);                  % Adds a reversible transport reaction for Ammonia.
model.rxns(end,1)={'TRP_c->e_Ammonia'};                                                     % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_Ammonia'))),:)=[5];           % Assigns a confidence score of 5.
model=changeGeneAssoc(model,'TRP_c->e_Ammonia','HPF_19190');  
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'in',lastMet);                                     % Adds a consuming exchange reaction for Ammonia.
model.rxns(end,1)={'EXC_IN_Ammonia'};                                                       % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_IN_Ammonia'))),:)=[5];             % Assigns a confidence score of 5.                                                                                                                                                   

                                                                            % Furthermore, producing exchange reactions are added for CO and CO2:

[model, addedRxns]=addExchangeRxns(model,'out',lastMetCO);                                  % Adds a producing exchange reaction for CO.
model.rxns(end,1)={'EXC_OUT_CO'};                                                           % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_OUT_CO'))),:)=[5];                 % Assigns a confidence score of 5.
[model, addedRxns]=addExchangeRxns(model,'out',lastMetCO2);                                 % Adds a producing exchange reaction for CO2. 
model.rxns(end,1)={'EXC_OUT_CO2'};                                                          % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_OUT_CO2'))),:)=[5];                % Assigns a confidence score of 5.

                                                                            % Finally, although not experimentally verified, producing exchange reactions were added for 
                                                                            % Acetate and Glycerol as excretion of these compounds are commonly observed in other 
                                                                            % microorganisms and often occur as a consequence of cells having to re-oxidize NADH and NADPH:                                                                            

[model, addedRxns]=addExchangeRxns(model,'out',lastMetAcetate);                             % Adds a producing exchange reaction for Acetate. 
model.rxns(end,1)={'EXC_OUT_Acetate'};                                                      % Assigns an appropriate reaction ID to the exchange reaction.                                                                        
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_OUT_Acetate'))),:)=[2];            % Assigns a confidence score of 2 (due to lack of evidence) to the producing exchange
                                                                                            % reaction for Acetate.
[model, addedRxns]=addExchangeRxns(model,'out',lastMetGlycerol);                            % Adds a producing exchange reaction for Glycerol.  
model.rxns(end,1)={'EXC_OUT_Glycerol'};                                                     % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_OUT_Glycerol'))),:)=[2];           % Assigns a confidence score of 2 (due to lack of evidence) to the producing exchange
                                                                                            % reaction for Glycerol.


%% Manual curation of reaction directionalities pertaining to reactions automatically generated in the FDR
                                                   
                                                                            % All reaction directionalities are - unless otherwise stated - curated using the MetaCyc database
                                                                            % (Caspi et al., 2017) where directionalities have been curated manually by trained scientists. URLs
                                                                            % for each reaction point to specific entries in the MetaCyc database from which information on 
                                                                            % reaction directionality was taken (when available).
                                                                            
                                                                            % Note also that the reactants of some reactions have been swapped to products and vice versa
                                                                            % through commands utilizing the changeRxns-function below. This is done so as to be able to
                                                                            % properly constrain the reaction bounds.
                                                
                                                                            % Reaction directionalities are set by manipulating the reaction bounds accordingly:                                                                                                                        
                                                                            
model=setParam(model,'lb',{'R00005'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ALLOPHANATE-HYDROLASE-RXN&redirect=T
model=changeRxns(model,'R00006','2 Pyruvate[c] <=> CO2[c] + 2-Acetolactate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00006'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETOLACTSYN-RXN
%R00008 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=4.1.3.17-RXN&redirect=T
model=setParam(model,'lb',{'R00013'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLYOCARBOLIG-RXN&redirect=T
model=setParam(model,'lb',{'R00014'},[0]);model=setParam(model,'ub',{'R00014'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12583&redirect=T
model=setParam(model,'lb',{'R00022'},[0]);                                                  %https://biocyc.org/GCF_004353865/NEW-IMAGE?type=REACTION&object=RXN-12625
model=setParam(model,'lb',{'R00024'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN&redirect=T
model=setParam(model,'lb',{'R00025'},[0]);model=setParam(model,'ub',{'R00025'},[1000]);     %https://biocyc.org/GCF_004353865/NEW-IMAGE?type=REACTION&object=2-NITROPROPANE-DIOXYGENASE-RXN
model=setParam(model,'lb',{'R00026'},[0]);                                                  %https://biocyc.org/ECOLI/NEW-IMAGE?type=REACTION&object=3.2.1.21-RXN
model=setParam(model,'lb',{'R00028'},[0]);                                                  %https://biocyc.org/ECOLI/NEW-IMAGE?type=REACTION&object=MALTODEXGLUCOSID-RXN
model=setParam(model,'lb',{'R00036'},[0]);                                                  %https://biocyc.org/ECOLI/NEW-IMAGE?type=REACTION&object=PORPHOBILSYNTH-RXN
model=setParam(model,'lb',{'R00066'},[0]);                                                  %https://biocyc.org/ECOLI/NEW-IMAGE?type=REACTION&object=RIBOFLAVIN-SYN-RXN
model=setParam(model,'lb',{'R00078'},[0]);model=setParam(model,'ub',{'R00078'},[1000]);     %https://biocyc.org/ECOLI/NEW-IMAGE?type=REACTION&object=RXN0-1483
model=setParam(model,'lb',{'R00081'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=CYTOCHROME-C-OXIDASE-RXN&redirect=T
model=setParam(model,'lb',{'R00084'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=OHMETHYLBILANESYN-RXN
model=setParam(model,'lb',{'R00089'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENYLATECYC-RXN
model=changeRxns(model,'R00094','NADH[c] + H+[c] + Glutathione disulfide[c] <=> NAD+[c] + 2 Glutathione[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00094'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTATHIONE-REDUCT-NADPH-RXN
model=setParam(model,'lb',{'R00104'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NAD-KIN-RXN
%R00112 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PYRNUTRANSHYDROGEN-RXN
model=changeRxns(model,'R00114','NADPH[c] + 2-Oxoglutarate[c] + L-Glutamine[c] + H+[c] <=> NADP+[c] + 2 L-Glutamate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00114'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTAMATESYN-RXN
model=changeRxns(model,'R00115','NADPH[c] + H+[c] + Glutathione disulfide[c] <=> NADP+[c] + 2 Glutathione[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00115'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTATHIONE-REDUCT-NADPH-RXN
model=setParam(model,'lb',{'R00125'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3.6.1.41-RXN
%R00127 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENYL-KIN-RXN
model=setParam(model,'lb',{'R00130'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DEPHOSPHOCOAKIN-RXN
%R00131 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UREASE-RXN
%R00132 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CARBODEHYDRAT-RXN
model=setParam(model,'lb',{'R00137'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2.7.7.1-RXN
model=setParam(model,'lb',{'R00139'},[0]);model=setParam(model,'ub',{'R00139'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NUCLEOSIDE-DIP-KIN-RXN
model=changeRxns(model,'R00143','NADH[c] + H+[c] + Hydroxylamine[c] <=> H2O[c] + NAD+[c] + Ammonia[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00143'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HYDROXYLAMINE-REDUCTASE-NADH-RXN
model=setParam(model,'lb',{'R00148'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=AMONITRO-RXN
model=setParam(model,'lb',{'R00156'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UDPKIN-RXN
model=setParam(model,'lb',{'R00158'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12002
model=setParam(model,'lb',{'R00161'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=FADSYN-RXN
model=changeRxns(model,'R00177','H2O[c] + ATP[c] + L-Methionine[c] <=> Orthophosphate[c] + Diphosphate[c] + S-Adenosyl-L-methionine[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00177'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=S-ADENMETSYN-RXN
model=setParam(model,'lb',{'R00178'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SAMDECARB-RXN
model=setParam(model,'lb',{'R00182'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=AMP-NUCLEOSID-RXN
model=setParam(model,'lb',{'R00183'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=AMP-DEPHOSPHORYLATION-RXN
model=setParam(model,'lb',{'R00185'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENOSINE-KINASE-RXN
model=setParam(model,'lb',{'R00188'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=325-BISPHOSPHATE-NUCLEOTIDASE-RXN
model=changeRxns(model,'R00190','5-Phospho-alpha-D-ribose 1-diphosphate[c] + Adenine[c] <=> Diphosphate[c] + AMP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00190'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENPRIBOSYLTRAN-RXN
model=setParam(model,'lb',{'R00192'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENOSYLHOMOCYSTEINASE-RXN
model=setParam(model,'lb',{'R00194'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENOSYLHOMOCYSTEINE-NUCLEOSIDASE-RXN
model=setParam(model,'lb',{'R00196'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=L-LACTATE-DEHYDROGENASE-CYTOCHROME-RXN
model=setParam(model,'lb',{'R00197'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=D-LACTATE-DEHYDROGENASE-CYTOCHROME-RXN
%R00199 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PEPSYNTH-RXN
model=changeRxns(model,'R00200','ADP[c] + Phosphoenolpyruvate[c] <=> ATP[c] + Pyruvate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00200'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PEPDEPHOS-RXN
model=setParam(model,'lb',{'R00209'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-1133
model=setParam(model,'lb',{'R00215'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.1.1.83-RXN
model=setParam(model,'lb',{'R00216'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MALIC-NADP-RXN
model=setParam(model,'lb',{'R00220'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4.3.1.17-RXN
model=setParam(model,'lb',{'R00221'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DSERDEAM-RXN
model=changeRxns(model,'R00226','2 Pyruvate[c] <=> CO2[c] + (S)-2-Acetolactate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00226'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETOLACTSYN-RXN
%R00228 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETALD-DEHYDROG-RXN
model=setParam(model,'lb',{'R00233'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MALONYL-COA-DECARBOXYLASE-RXN
model=setParam(model,'lb',{'R00235'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETATE--COA-LIGASE-RXN
model=setParam(model,'lb',{'R00236'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETATE--COA-LIGASE-RXN
%R00238 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETYL-COA-ACETYLTRANSFER-RXN
model=setParam(model,'lb',{'R00239'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTKIN-RXN
model=setParam(model,'lb',{'R00243'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTAMATE-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R00245'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14116
%R00248 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTAMATE-DEHYDROGENASE-NADP%2b-RXN
model=setParam(model,'lb',{'R00251'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=5-OXOPROLINASE-ATP-HYDROLYSING-RXN
model=setParam(model,'lb',{'R00253'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTAMINESYN-RXN
model=setParam(model,'lb',{'R00256'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTAMIN-RXN
model=setParam(model,'lb',{'R00257'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NAD-SYNTH-GLN-RXN
%R00258 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ALANINE-AMINOTRANSFERASE-RXN
%R00259 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=N-ACETYLTRANSFER-RXN
%R00260 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTRACE-RXN
model=setParam(model,'lb',{'R00264'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=25-DIOXOVALERATE-DEHYDROGENASE-RXN
%R00267 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ISOCITDEH-RXN
%R00268 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8642
model=setParam(model,'lb',{'R00274'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTATHIONE-PEROXIDASE-RXN
model=setParam(model,'lb',{'R00277'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PMPOXI-RXN
model=setParam(model,'lb',{'R00278'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PNPOXI-RXN
model=setParam(model,'lb',{'R00286'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UGD-RXN
%R00289 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUC1PURIDYLTRANS-RXN
%R00291 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UDPGLUCEPIM-RXN
%R00294 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=NITRIC-OXIDE-REDUCTASE-RXN&redirect=T
model=setParam(model,'lb',{'R00299'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUCOKIN-RXN
model=setParam(model,'lb',{'R00306'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-10773
%R00310 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PROTOHEMEFERROCHELAT-RXN
model=setParam(model,'lb',{'R00316'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETATE--COA-LIGASE-RXN
model=setParam(model,'lb',{'R00330'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GDPKIN-RXN
%R00332 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GUANYL-KIN-RXN
model=setParam(model,'lb',{'R00336'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PPGPPSYN-RXN
%R00342 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MALATE-DEH-RXN
model=changeRxns(model,'R00345','H2O[c] + CO2[c] + Phosphoenolpyruvate[c] <=> Orthophosphate[c] + Oxaloacetate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00345'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PEPCARBOX-RXN
%R00350 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-2464
model=changeRxns(model,'R00351','H2O[c] + Acetyl-CoA[c] + Oxaloacetate[c] <=> CoA[c] + Citrate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00351'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CITSYN-RXN
model=setParam(model,'lb',{'R00357'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=L-AMINO-ACID-OXIDASE-RXN
model=setParam(model,'lb',{'R00362'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CITLY-RXN
%R00401 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ALARACECAT-RXN
%R00405 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SUCCCOASYN-RXN
%R00410 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXNI-2
model=setParam(model,'lb',{'R00416'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NAG1P-URIDYLTRANS-RXN
model=setParam(model,'lb',{'R00425'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GTP-CYCLOHYDRO-II-RXN
model=setParam(model,'lb',{'R00426'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14140
                                    model=setParam(model,'lb',{'R00428'},[-1000]);model=setParam(model,'ub',{'R00428'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00428'))),:)=[3];  % Reaction directionality unclear!                                                                     
model=setParam(model,'lb',{'R00429'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GTPPYPHOSKIN-RXN
model=changeRxns(model,'R00430','GDP[c] + Phosphoenolpyruvate[c] <=> Pyruvate[c] + GTP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00430'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14117
%R00431 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4.1.1.32-RXN
model=setParam(model,'lb',{'R00434'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GUANYLCYC-RXN
model=setParam(model,'lb',{'R00449'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=LYSINE-2-MONOOXYGENASE-RXN
model=setParam(model,'lb',{'R00451'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DIAMINOPIMDECARB-RXN
model=changeRxns(model,'R00465','NADPH[c] + Glyoxylate[c] + H+[c] <=> NADP+[c] + Glycolate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00465'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLYOXYLATE-REDUCTASE-NADP%2b-RXN
%R00470 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4OH2OXOGLUTARALDOL-RXN
model=changeRxns(model,'R00472','H2O[c] + Acetyl-CoA[c] + Glyoxylate[c] <=> CoA[c] + (S)-Malate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00472'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MALSYN-RXN
model=setParam(model,'lb',{'R00475'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-969
%R00479 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ISOCIT-CLEAV-RXN
model=setParam(model,'lb',{'R00480'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ASPARTATEKIN-RXN
model=setParam(model,'lb',{'R00481'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=L-ASPARTATE-OXID-RXN
model=setParam(model,'lb',{'R00485'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ASPARAGHYD-RXN
model=setParam(model,'lb',{'R00489'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ASPDECARBOX-RXN
model=setParam(model,'lb',{'R00491'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ASPARTATE-RACEMASE-RXN
model=setParam(model,'lb',{'R00494'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12618
model=setParam(model,'lb',{'R00497'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTATHIONE-SYN-RXN
model=setParam(model,'lb',{'R00511'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14026
%R00519 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.2.1.2-RXN
model=setParam(model,'lb',{'R00527'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=S-FORMYLGLUTATHIONE-HYDROLASE-RXN
%R00529 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SULFATE-ADENYLYLTRANS-RXN
model=setParam(model,'lb',{'R00548'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-5187
model=setParam(model,'lb',{'R00549'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RIBOFLAVINKIN-RXN
%R00551 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ARGINASE-RXN
%R00566 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ARGDECARBOX-RXN
model=setParam(model,'lb',{'R00566'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14118
model=setParam(model,'lb',{'R00570'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CDPKIN-RXN
model=setParam(model,'lb',{'R00571'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14325
model=setParam(model,'lb',{'R00573'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CTPSYN-RXN
model=setParam(model,'lb',{'R00575'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CARBPSYN-RXN
model=setParam(model,'lb',{'R00578'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ASNSYNB-RXN
model=setParam(model,'lb',{'R00582'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-5114
model=setParam(model,'lb',{'R00586'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SERINE-O-ACETTRAN-RXN
model=setParam(model,'lb',{'R00597'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12458
model=setParam(model,'lb',{'R00601'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-5144
model=setParam(model,'lb',{'R00602'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14189
model=setParam(model,'lb',{'R00606'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-2841
%R00615 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-3542
model=setParam(model,'lb',{'R00617'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THI-P-KIN-RXN
                                    model=setParam(model,'lb',{'R00621'},[-1000]);model=setParam(model,'ub',{'R00621'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00621'))),:)=[3];  % Reaction directionality unclear!
%R00658 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2PGADEHYDRAT-RXN
model=setParam(model,'lb',{'R00660'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UDPNACETYLGLUCOSAMENOLPYRTRANS-RXN
model=setParam(model,'lb',{'R00674'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-2382
model=setParam(model,'lb',{'R00678'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8665
model=setParam(model,'lb',{'R00691'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CARBOXYCYCLOHEXADIENYL-DEHYDRATASE-RXN
%R00694 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-10814
%R00700 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HYDROGEN-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R00705'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-2902
model=setParam(model,'lb',{'R00706'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9958
model=setParam(model,'lb',{'R00707'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PYRROLINECARBDEHYDROG-RXN
model=setParam(model,'lb',{'R00708'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PYRROLINECARBDEHYDROG-RXN
model=setParam(model,'lb',{'R00710'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN66-3
model=setParam(model,'lb',{'R00711'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-3962
model=setParam(model,'lb',{'R00713'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SUCCINATE-SEMIALDEHYDE-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R00714'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SUCCSEMIALDDEHYDROG-RXN
%R00717 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLYCOLATE-REDUCTASE-RXN
model=setParam(model,'lb',{'R00720'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-6382
%R00722 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14120
%R00726 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12481
%R00734 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=TYROSINE-AMINOTRANSFERASE-RXN
model=setParam(model,'lb',{'R00742'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETYL-COA-CARBOXYLTRANSFER-RXN
model=setParam(model,'lb',{'R00749'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ETHAMLY-RXN
%R00750 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MHPELY-RXN
%R00751 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THREONINE-ALDOLASE-RXN
%R00754 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ALCOHOL-DEHYDROG-RXN
model=setParam(model,'lb',{'R00756'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=6PFRUCTPHOS-RXN
model=setParam(model,'lb',{'R00760'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=FRUCTOKINASE-RXN
model=setParam(model,'lb',{'R00762'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=F16BDEPHOS-RXN
%R00764 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2.7.1.90-RXN
%R00768 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=L-GLN-FRUCT-6-P-AMINOTRANS-RXN
%R00771 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PGLUCISOM-RXN
model=setParam(model,'lb',{'R00774'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UREA-CARBOXYLASE-RXN
model=setParam(model,'lb',{'R00782'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=LCYSDESULF-RXN
model=changeRxns(model,'R00783','H+[c] + Nitrite[c] + Ferrocytochrome c[c] <=> H2O[c] + Ferricytochrome c[c] + Nitric oxide[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00783'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NITRITE-REDUCTASE-CYTOCHROME-RXN
model=changeRxns(model,'R00785','H+[c] + Nitrite[c] + Reduced azurin[c] <=> H2O[c] + Nitric oxide[c] + Oxidized azurin[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00785'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NITRITE-REDUCTASE-CYTOCHROME-RXN
model=changeRxns(model,'R00787','3 NADH[c] + 3 H+[c] + Nitrite[c] <=> 2 H2O[c] + 3 NAD+[c] + Ammonia[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00787'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13854
                                    model=setParam(model,'lb',{'R00798'},[-1000]);model=setParam(model,'ub',{'R00798'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00798'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R00801'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3.2.1.48-RXN
model=setParam(model,'lb',{'R00802'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3.2.1.48-RXN
model=setParam(model,'lb',{'R00813'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.3.1.25-RXN    
model=setParam(model,'lb',{'R00815'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHENOL-2-MONOOXYGENASE-RXN
model=setParam(model,'lb',{'R00816'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CATECHOL-23-DIOXYGENASE-RXN
model=setParam(model,'lb',{'R00818'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SALICYLATE-1-MONOOXYGENASE-RXN
model=changeRxns(model,'R00829','CoA[c] + 3-Oxoadipyl-CoA[c] <=> Acetyl-CoA[c] + Succinyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00829'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-3641
%R00833 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=METHYLMALONYL-COA-MUT-RXN
model=setParam(model,'lb',{'R00835'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLU6PDEHYDROG-RXN
model=changeRxns(model,'R00842','NADH[c] + H+[c] + Glycerone phosphate[c] <=> NAD+[c] + sn-Glycerol 3-phosphate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00842'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.1.1.8-RXN
model=changeRxns(model,'R00844','NADPH[c] + H+[c] + Glycerone phosphate[c] <=> NADP+[c] + sn-Glycerol 3-phosphate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00844'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLYC3PDEHYDROGBIOSYN-RXN
%R00847 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLYCEROL-KIN-RXN
model=setParam(model,'lb',{'R00848'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.1.5.3&quickSearch=Quick+Search
model=changeRxns(model,'R00858','3 NADPH[c] + 3 H+[c] + Sulfite[c] <=> 3 H2O[c] + 3 NADP+[c] + Hydrogen sulfide[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00858'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SULFITE-REDUCT-RXN
model=setParam(model,'lb',{'R00867'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=FRUCTOKINASE-RXN
model=setParam(model,'lb',{'R00868'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MANNITOL-2-DEHYDROGENASE-RXN
%R00878 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUCISOM-RXN
model=setParam(model,'lb',{'R00885'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2.7.7.13-RXN
model=setParam(model,'lb',{'R00888'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GDPMANDEHYDRA-RXN
model=setParam(model,'lb',{'R00894'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTCYSLIG-RXN
%R00897 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACSERLY-RXN
model=setParam(model,'lb',{'R00899'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-6622
model=setParam(model,'lb',{'R00904'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-6382
model=setParam(model,'lb',{'R00905'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=BETA-UREIDOPROPIONASE-RXN
%R00907 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2.6.1.18-RXN
%R00908 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-2901
model=setParam(model,'lb',{'R00922'},[0]);model=setParam(model,'ub',{'R00922'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.2.1.27-RXN
model=setParam(model,'lb',{'R00925'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PROPIONATE--COA-LIGASE-RXN
                                    model=setParam(model,'lb',{'R00926'},[-1000]);model=setParam(model,'ub',{'R00926'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00926'))),:)=[3];  % Reaction directionality unclear!
%R00927 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=METHYLACETOACETYLCOATHIOL-RXN
model=setParam(model,'lb',{'R00935'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11213
model=changeRxns(model,'R00936','NADH[c] + H+[c] + Dihydrofolate[c] <=> NAD+[c] + Tetrahydrofolate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00936'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DIHYDROFOLATEREDUCT-RXN
model=changeRxns(model,'R00937','2 NADH[c] + 2 H+[c] + Folate[c] <=> 2 NAD+[c] + Tetrahydrofolate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00937'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.5.1.3&quickSearch=Quick+Search
model=changeRxns(model,'R00939','NADPH[c] + H+[c] + Dihydrofolate[c] <=> NADP+[c] + Tetrahydrofolate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00939'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.5.1.3&quickSearch=Quick+Search
model=changeRxns(model,'R00940','2 NADPH[c] + 2 H+[c] + Folate[c] <=> 2 NADP+[c] + Tetrahydrofolate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00940'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.5.1.3&quickSearch=Quick+Search
model=setParam(model,'lb',{'R00942'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-6341
model=setParam(model,'lb',{'R00944'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=FORMYLTHFDEFORMYL-RXN
%R00945 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLYOHMETRANS-RXN
model=setParam(model,'lb',{'R00946'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HOMOCYSMETB12-RXN
model=setParam(model,'lb',{'R00948'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUC1PADENYLTRANS-RXN
%R00956 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2.7.7.33-RXN
%R00959 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHOSPHOGLUCMUT-RXN
model=setParam(model,'lb',{'R00963'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14025
model=setParam(model,'lb',{'R00965'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=OROTPDECARB-RXN
model=changeRxns(model,'R00966','Uracil[c] + 5-Phospho-alpha-D-ribose 1-diphosphate[c] <=> Diphosphate[c] + UMP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.
model=setParam(model,'lb',{'R00966'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=URACIL-PRIBOSYLTRANS-RXN
model=setParam(model,'lb',{'R00974'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CYTDEAM-RXN
model=setParam(model,'lb',{'R00982'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=AMINOBENZCOALIG-RXN
                                    model=setParam(model,'lb',{'R00985'},[-1000]);model=setParam(model,'ub',{'R00985'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00985'))),:)=[3];  % Reaction directionality unclear!
%R00986 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ANTHRANSYN-RXN
model=setParam(model,'lb',{'R00987'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=KYNURENINASE-RXN
                                    model=setParam(model,'lb',{'R00988'},[-1000]);model=setParam(model,'ub',{'R00988'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00988'))),:)=[3];  % Reaction directionality unclear!
model=changeRxns(model,'R00994','NAD+[c] + D-erythro-3-Methylmalate[c] <=> NADH[c] + CO2[c] + H+[c] + 2-Oxobutanoate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R00994'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.85&quickSearch=Quick+Search
model=setParam(model,'lb',{'R00996'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THREDEHYD-RXN
%R01015 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=TRIOSEPISOMERIZATION-RXN
model=setParam(model,'lb',{'R01016'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=METHGLYSYN-RXN
model=setParam(model,'lb',{'R01025'},[-1000]);model=setParam(model,'ub',{'R01025'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=CHD-RXN&redirect=T                                       Irreversible in the original FDR.
model=setParam(model,'lb',{'R01030'},[0]);model=setParam(model,'ub',{'R01030'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3.1.4.2-RXN
%R01049 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PRPPSYN-RXN
model=setParam(model,'lb',{'R01054'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-1441
%R01056 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RIB5PISOM-RXN
%R01057 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PPENTOMUT-RXN
%R01061 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GAPOXNPHOSPHN-RXN
%R01064 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DEHYDDEOXPHOSGALACT-ALDOL-RXN
%R01067 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2TRANSKETO-RXN
%R01068 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=F16ALDOLASE-RXN
%R01070 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=F16ALDOLASE-RXN
%R01071 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ATPPHOSPHORIBOSYLTRANS-RXN
%R01072 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PRPPAMIDOTRANS-RXN
%R01073 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PRTRANS-RXN
model=setParam(model,'lb',{'R01074'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-7192
model=setParam(model,'lb',{'R01078'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2.8.1.6-RXN
%R01082 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=FUMHYDR-RXN
%R01083 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=AMPSYN-RXN
model=setParam(model,'lb',{'R01085'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-10445
%R01086 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ARGSUCCINLYA-RXN
%R01090 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=BRANCHED-CHAINAMINOTRANSFERLEU-RXN
                                    model=setParam(model,'lb',{'R01106'},[-1000]);model=setParam(model,'ub',{'R01106'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R01106'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R01122'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-6274
model=setParam(model,'lb',{'R01126'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-7607
%R01127 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=IMPCYCLOHYDROLASE-RXN
%R01130 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=IMP-DEHYDROG-RXN
model=changeRxns(model,'R01134','NADPH[c] + H+[c] + GMP[c] <=> NADP+[c] + Ammonia[c] + IMP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01134'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GMP-REDUCT-RXN
model=setParam(model,'lb',{'R01135'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENYLOSUCCINATE-SYNTHASE-RXN
model=setParam(model,'lb',{'R01137'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DADPKIN-RXN
model=changeRxns(model,'R01138','Phosphoenolpyruvate[c] + dADP[c] <=> Pyruvate[c] + dATP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01138'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14192
model=setParam(model,'lb',{'R01150'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DALADALALIG-RXN
model=setParam(model,'lb',{'R01157'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=AGMATIN-RXN
model=setParam(model,'lb',{'R01158'},[0]);model=setParam(model,'ub',{'R01158'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8001
model=setParam(model,'lb',{'R01163'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HISTALDEHYD-RXN
model=setParam(model,'lb',{'R01168'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HISTIDINE-AMMONIA-LYASE-RXN
%R01172 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=BUTANAL-DEHYDROGENASE-RXN  
%R01177 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12565
model=setParam(model,'lb',{'R01185'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-5408
model=setParam(model,'lb',{'R01186'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-10952
model=setParam(model,'lb',{'R01187'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MYO-INOSITOL-1OR-4-MONOPHOSPHATASE-RXN
model=setParam(model,'lb',{'R01209'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DIHYDROXYISOVALDEHYDRAT-RXN
model=changeRxns(model,'R01213','H2O[c] + Acetyl-CoA[c] + 3-Methyl-2-oxobutanoic acid[c] <=> CoA[c] + alpha-Isopropylmalate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01213'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2-ISOPROPYLMALATESYN-RXN
%R01214 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=BRANCHED-CHAINAMINOTRANSFERVAL-RXN
%R01215 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=VALINE-PYRUVATE-AMINOTRANSFER-RXN
%R01220 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=METHYLENETHFDEHYDROG-NADP-RXN
%R01221 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GCVMULTI-RXN
model=changeRxns(model,'R01224','NADPH[c] + H+[c] + 5,10-Methylenetetrahydrofolate[c] <=> NADP+[c] + 5-Methyltetrahydrofolate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01224'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.5.1.20-RXN
%R01226 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3-CH3-2-OXOBUTANOATE-OH-CH3-XFER-RXN
model=setParam(model,'lb',{'R01227'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-7609
model=changeRxns(model,'R01229','5-Phospho-alpha-D-ribose 1-diphosphate[c] + Guanine[c] <=> Diphosphate[c] + GMP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01229'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GUANPRIBOSYLTRAN-RXN
model=setParam(model,'lb',{'R01230'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GMP-SYN-NH3-RXN
model=setParam(model,'lb',{'R01231'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GMP-SYN-GLUT-RXN
model=setParam(model,'lb',{'R01238'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4-HYDROXYBENZOATE-DECARBOXYLASE-RXN
model=setParam(model,'lb',{'R01244'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENINE-DEAMINASE-RXN
model=changeRxns(model,'R01248','NADH[c] + H+[c] + (S)-1-Pyrroline-5-carboxylate[c] <=> NAD+[c] + L-Proline[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01248'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PYRROLINECARBREDUCT-RXN
model=changeRxns(model,'R01251','NADPH[c] + H+[c] + (S)-1-Pyrroline-5-carboxylate[c] <=> NADP+[c] + L-Proline[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01251'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PYRROLINECARBREDUCT-RXN
model=setParam(model,'lb',{'R01252'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN490-3641
model=setParam(model,'lb',{'R01253'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14903
model=setParam(model,'lb',{'R01262'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-6601
model=setParam(model,'lb',{'R01274'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.1.2.-&quickSearch=Quick+Search
model=setParam(model,'lb',{'R01280'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=6.2.1.3&quickSearch=Quick+Search
model=setParam(model,'lb',{'R01286'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CYSTATHIONINE-BETA-LYASE-RXN
%R01287 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETYLHOMOSER-CYS-RXN
%R01288 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9384   
model=setParam(model,'lb',{'R01298'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4-HYDROXYBENZOATE-3-MONOOXYGENASE-RXN
model=setParam(model,'lb',{'R01301'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3.1.2.23-RXN
%R01324 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14047
%R01325 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACONITATEDEHYDR-RXN
model=setParam(model,'lb',{'R01334'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GPH-RXN
model=setParam(model,'lb',{'R01354'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PROPIONATE--COA-LIGASE-RXN
model=setParam(model,'lb',{'R01357'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETOACETATE--COA-LIGASE-RXN
model=setParam(model,'lb',{'R01360'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HYDROXYMETHYLGLUTARYL-COA-LYASE-RXN
%R01361 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3-HYDROXYBUTYRATE-DEHYDROGENASE-RXN
model=changeRxns(model,'R01364','H2O[c] + 4-Fumarylacetoacetate[c] <=> Fumarate[c] + Acetoacetate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01364'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=FUMARYLACETOACETASE-RXN
model=setParam(model,'lb',{'R01372'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-10815
%R01373 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PREPHENATEDEHYDRAT-RXN
model=setParam(model,'lb',{'R01374'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11193
model=setParam(model,'lb',{'R01384'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UDP-GLUCURONATE-DECARBOXYLASE-RXN
%R01385 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UDP-GLUCURONATE-4-EPIMERASE-RXN
%R01388 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLYCERATE-DEHYDROGENASE-RXN
%R01392 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HYDROXYPYRUVATE-REDUCTASE-RXN
%R01394 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-305
%R01395 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14196&redirect=T
%R01397 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ASPCARBTRANS-RXN
%R01398 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ORNCARBAMTRANSFER-RXN
model=setParam(model,'lb',{'R01401'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=METHYLTHIOADENOSINE-NUCLEOSIDASE-RXN
model=changeRxns(model,'R01403','NADH[c] + H+[c] + trans-2,3-Dehydroacyl-[acyl-carrier protein][c] <=> NAD+[c] + Acyl-[acyl-carrier protein][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01403'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.9&quickSearch=Quick+Search
model=changeRxns(model,'R01404','NADPH[c] + H+[c] + trans-2,3-Dehydroacyl-[acyl-carrier protein][c] <=> NADP+[c] + Acyl-[acyl-carrier protein][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01404'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ENOYL-ACP-REDUCT-NADPH-RXN
model=setParam(model,'lb',{'R01411'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14197
model=setParam(model,'lb',{'R01424'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HIPPURATE-HYDROLASE-RXN
model=setParam(model,'lb',{'R01429'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=D-XYLOSE-1-DEHYDROGENASE-RXN
%R01432 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=XYLISOM-RXN
model=setParam(model,'lb',{'R01466'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THRESYN-RXN
model=setParam(model,'lb',{'R01470'},[0]);model=setParam(model,'ub',{'R01470'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14160
model=setParam(model,'lb',{'R01492'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.17&quickSearch=Quick+Search
model=setParam(model,'lb',{'R01507'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.13.11.14-RXN
%R01512 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHOSGLYPHOS-RXN
%R01513 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PGLYCDEHYDROG-RXN
%R01518 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15513
model=setParam(model,'lb',{'R01519'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUCONOLACT-RXN
%R01523 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHOSPHORIBULOKINASE-RXN
model=changeRxns(model,'R01525','NADPH[c] + H+[c] + D-Ribulose 5-phosphate[c] <=> NADP+[c] + D-Ribitol 5-phosphate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01525'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RIBITOL-5-PHOSPHATE-2-DEHYDROGENASE-RXN
%R01529 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RIBULP3EPIM-RXN
%R01530 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DARAB5PISOM-RXN
model=setParam(model,'lb',{'R01547'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DEOXYADENYLATE-KINASE-RXN
%R01561 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENPHOSPHOR-RXN
model=setParam(model,'lb',{'R01569'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THYMIDYLATE-5-PHOSPHATASE-RXN
%R01570 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THYM-PHOSPH-RXN
model=setParam(model,'lb',{'R01600'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUCOKIN-RXN
%R01602 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ALDOSE-1-EPIMERASE-RXN
model=setParam(model,'lb',{'R01625'},[0]);model=setParam(model,'ub',{'R01625'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HOLO-ACP-SYNTH-RXN
model=setParam(model,'lb',{'R01626'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MALONYL-COA-ACP-TRANSACYL-RXN
model=setParam(model,'lb',{'R01631'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PROTOCATECHUATE-34-DIOXYGENASE-RXN
model=setParam(model,'lb',{'R01632'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PROTOCATECHUATE-45-DIOXYGENASE-RXN
model=setParam(model,'lb',{'R01639'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=XYLULOKIN-RXN
%R01641 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1TRANSKETO-RXN
model=changeRxns(model,'R01644','NADH[c] + H+[c] + Succinate semialdehyde[c] <=> NAD+[c] + 4-Hydroxybutanoic acid[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01644'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4-HYDROXYBUTYRATE-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R01645'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4-HYDROXY-2-KETOPIMELATE-LYSIS-RXN
model=changeRxns(model,'R01647','2,4-Dihydroxyhept-2-enedioate[c] <=> Pyruvate[c] + Succinate semialdehyde[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01647'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4-HYDROXY-2-KETOPIMELATE-LYSIS-RXN
%R01648 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GABATRANSAM-RXN
%R01655 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=METHENYLTHFCYCLOHYDRO-RXN
model=setParam(model,'lb',{'R01658'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GPPSYN-RXN
model=setParam(model,'lb',{'R01664'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-5292
model=setParam(model,'lb',{'R01676'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GUANINE-DEAMINASE-RXN
model=setParam(model,'lb',{'R01687'},[0]);model=setParam(model,'ub',{'R01687'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GAMMA-GLUTAMYLTRANSFERASE-RXN
model=setParam(model,'lb',{'R01698'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.8.1.4-RXN&redirect=T
model=setParam(model,'lb',{'R01699'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-1134
model=setParam(model,'lb',{'R01710'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PMPOXI-RXN
model=setParam(model,'lb',{'R01711'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PYRIDOXINE-4-OXIDASE-RXN
model=setParam(model,'lb',{'R01714'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CHORISMATE-SYNTHASE-RXN
%R01715 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CHORISMATEMUT-RXN
model=setParam(model,'lb',{'R01716'},[-1000]);model=setParam(model,'ub',{'R01716'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PABASYN-RXN                                         Irreversible in the original FDR.
model=changeRxns(model,'R01724','H2O[c] + ATP[c] + H+[c] + 5-Phospho-alpha-D-ribose 1-diphosphate[c] + Nicotinate[c] <=> ADP[c] + Orthophosphate[c] + Diphosphate[c] + Nicotinate D-ribonucleotide[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01724'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NICOTINATEPRIBOSYLTRANS-RXN
model=setParam(model,'lb',{'R01728'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PREPHENATEDEHYDROG-RXN
%R01731 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PREPHENATE-ASP-TRANSAMINE-RXN
model=setParam(model,'lb',{'R01736'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLYOXII-RXN
model=setParam(model,'lb',{'R01737'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUCONOKIN-RXN
model=changeRxns(model,'R01739','NADPH[c] + H+[c] + 2-Keto-D-gluconic acid[c] <=> NADP+[c] + D-Gluconic acid[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01739'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.1.1.215-RXN
%R01751 is reversible (= in MetaCyc)                                                        %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TARTRATE-DECARBOXYLASE-RXN&redirect=T
model=setParam(model,'lb',{'R01752'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.2.1.3&quickSearch=Quick+Search
model=setParam(model,'lb',{'R01768'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-7682
model=setParam(model,'lb',{'R01771'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HOMOSERKIN-RXN
model=changeRxns(model,'R01773','NADH[c] + H+[c] + L-Aspartate 4-semialdehyde[c] <=> NAD+[c] + L-Homoserine[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01773'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HOMOSERDEHYDROG-RXN
model=changeRxns(model,'R01775','NADPH[c] + H+[c] + L-Aspartate 4-semialdehyde[c] <=> NADP+[c] + L-Homoserine[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01775'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HOMOSERDEHYDROG-RXN
%R01776 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HOMOSERINE-O-ACETYLTRANSFERASE-RXN
model=setParam(model,'lb',{'R01777'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=HOMSUCTRAN-RXN
model=changeRxns(model,'R01779','NADPH[c] + H+[c] + 3-Oxoacyl-CoA[c] <=> NADP+[c] + (3R)-3-Hydroxyacyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01779'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-7698&redirect=T
model=setParam(model,'lb',{'R01786'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUCOKIN-RXN
model=setParam(model,'lb',{'R01795'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.16.1&quickSearch=Quick+Search
model=setParam(model,'lb',{'R01799'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CDPDIGLYSYN-RXN
model=setParam(model,'lb',{'R01800'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHOSPHASERSYN-RXN
model=setParam(model,'lb',{'R01801'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHOSPHAGLYPSYN-RXN
%R01818 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHOSMANMUT-RXN
%R01819 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MANNPISOM-RXN
%R01826 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DAHPSYN-RXN
%R01827 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=TRANSALDOL-RXN
%R01829 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=4.1.2.13&quickSearch=Quick+Search
%R01830 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2TRANSKETO-RXN
model=setParam(model,'lb',{'R01845'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SEDOHEPTULOSE-BISPHOSPHATASE-RXN
model=setParam(model,'lb',{'R01855'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-385
model=setParam(model,'lb',{'R01856'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DGTPTRIPHYDRO-RXN
model=setParam(model,'lb',{'R01857'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DGDPKIN-RXN
model=changeRxns(model,'R01858','Phosphoenolpyruvate[c] + dGDP[c] <=> Pyruvate[c] + dGTP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01858'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14207   
%R01859 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PROPIONYL-COA-CARBOXY-RXN
%R01863 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=INOPHOSPHOR-RXN
model=setParam(model,'lb',{'R01868'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.5.2&quickSearch=Quick+Search
%R01870 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=OROPRIBTRANS-RXN
%R01876 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=URPHOS-RXN
%R01884 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CREATININASE-RXN
%R01895 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RIBITOL-2-DEHYDROGENASE-RXN
%R01899 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9951
%R01900 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACONITATEHYDR-RXN
%R01920 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SPERMIDINESYN-RXN
model=setParam(model,'lb',{'R01931'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THIOSULFATE-SULFURTRANSFERASE-RXN
model=setParam(model,'lb',{'R01933'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2-KETO-ADIPATE-DEHYDROG-RXN&redirect=T
%R01939 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2-AMINOADIPATE-AMINOTRANSFERASE-RXN
model=setParam(model,'lb',{'R01954'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ARGSUCCINSYN-RXN
model=setParam(model,'lb',{'R01959'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ARYLFORMAMIDASE-RXN
model=setParam(model,'lb',{'R01968'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14142
%R01975 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11662
model=changeRxns(model,'R01976','NADPH[c] + H+[c] + Acetoacetyl-CoA[c] <=> NADP+[c] + (S)-3-Hydroxybutanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01976'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3-HYDROXYBUTYRYL-COA-DEHYDROGENASE-RXN
model=changeRxns(model,'R01977','NADPH[c] + H+[c] + Acetoacetyl-CoA[c] <=> NADP+[c] + (R)-3-Hydroxybutanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R01977'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-5901
%R01986 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14209&redirect=T
model=setParam(model,'lb',{'R01990'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GUANIDINOBUTYRASE-RXN
%R01993 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DIHYDROOROT-RXN
model=setParam(model,'lb',{'R02003'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=FPPSYN-RXN
model=changeRxns(model,'R02016','NADPH[c] + H+[c] + Thioredoxin disulfide[c] <=> NADP+[c] + Thioredoxin[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02016'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THIOREDOXIN-REDUCT-NADPH-RXN
model=changeRxns(model,'R02017','ADP[c] + Thioredoxin[c] <=> H2O[c] + dADP[c] + Thioredoxin disulfide[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02017'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADPREDUCT-RXN
model=changeRxns(model,'R02018','UDP[c] + Thioredoxin[c] <=> H2O[c] + Thioredoxin disulfide[c] + dUDP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02018'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UDPREDUCT-RXN
model=changeRxns(model,'R02019','GDP[c] + Thioredoxin[c] <=> H2O[c] + Thioredoxin disulfide[c] + dGDP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02019'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GDPREDUCT-RXN
%R02021 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.8.4.8-RXN
model=changeRxns(model,'R02024','CDP[c] + Thioredoxin[c] <=> H2O[c] + Thioredoxin disulfide[c] + dCDP[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02024'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CDPREDUCT-RXN
model=setParam(model,'lb',{'R02026'},[0]);                                                  %https://www.genome.jp/dbget-bin/www_bget?rn:R02026
model=setParam(model,'lb',{'R02029'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.1.3.27&quickSearch=Quick+Search
model=setParam(model,'lb',{'R02035'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=6PGLUCONOLACT-RXN
model=setParam(model,'lb',{'R02036'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PGLUCONDEHYDRAT-RXN
model=setParam(model,'lb',{'R02055'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHOSPHASERDECARB-RXN
%R02060 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=5.4.2.10-RXN
model=setParam(model,'lb',{'R02065'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2.5.1.32-RXN
%R02073 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2.7.1.90-RXN
%R02085 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=METHYLGLUTACONYL-COA-HYDRATASE-RXN
model=setParam(model,'lb',{'R02088'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14161
%R02090 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GMKALT-RXN
model=setParam(model,'lb',{'R02093'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DTDPKIN-RXN&redirect=T
model=setParam(model,'lb',{'R02094'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DTMPKI-RXN
model=setParam(model,'lb',{'R02098'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14122
model=setParam(model,'lb',{'R02100'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DUTP-PYROP-RXN
model=setParam(model,'lb',{'R02101'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THYMIDYLATESYN-RXN
model=setParam(model,'lb',{'R02102'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14143
%R02103 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-901
model=setParam(model,'lb',{'R02110'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-7710&redirect=T
model=setParam(model,'lb',{'R02111'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.4.1.1&quickSearch=Quick+Search
%R02124 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RETINOL-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R02135'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXNQT-4191
%R02147 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-5199
model=setParam(model,'lb',{'R02164'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.5.1&quickSearch=Quick+Search
%R02199 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=BRANCHED-CHAINAMINOTRANSFERILEU-RXN
model=setParam(model,'lb',{'R02222'},[0]);model=setParam(model,'ub',{'R02222'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.14.19.1-RXN
model=changeRxns(model,'R02235','NADH[c] + H+[c] + Folate[c] <=> NAD+[c] + Dihydrofolate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02235'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18357
model=changeRxns(model,'R02236','NADPH[c] + H+[c] + Folate[c] <=> NADP+[c] + Dihydrofolate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02236'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18357
model=setParam(model,'lb',{'R02237'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DIHYDROFOLATESYNTH-RXN
model=setParam(model,'lb',{'R02240'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DIACYLGLYKIN-RXN
model=changeRxns(model,'R02241','Acyl-CoA[c] + 1-Acyl-sn-glycerol 3-phosphate[c] <=> CoA[c] + Phosphatidate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02241'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-1623
model=changeRxns(model,'R02251','Acyl-CoA[c] + 1,2-Diacyl-sn-glycerol[c] <=> CoA[c] + Triacylglycerol[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02251'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DIACYLGLYCEROL-O-ACYLTRANSFERASE-RXN
model=changeRxns(model,'R02272','(S)-4-Amino-5-oxopentanoate[c] <=> 5-Aminolevulinate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02272'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GSAAMINOTRANS-RXN
model=setParam(model,'lb',{'R02273'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=5-AMINOPENTANAMIDASE-RXN
%R02274 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=VAGL-RXN
model=setParam(model,'lb',{'R02278'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4.2.1.43-RXN
model=setParam(model,'lb',{'R02279'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=4.2.1.41-RXN
model=setParam(model,'lb',{'R02282'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTAMATE-N-ACETYLTRANSFERASE-RXN
%R02283 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETYLORNTRANSAM-RXN
model=setParam(model,'lb',{'R02285'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=FORMIMINOGLUTAMASE-RXN
model=setParam(model,'lb',{'R02288'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=IMIDAZOLONEPROPIONASE-RXN
%R02291 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ASPARTATE-SEMIALDEHYDE-DEHYDROGENASE-RXN
%R02296 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=2.4.2.2&quickSearch=Quick+Search
model=setParam(model,'lb',{'R02297'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=XANTHOSINEPHOSPHORY-RXN
model=setParam(model,'lb',{'R02300'},[0]);model=setParam(model,'ub',{'R02300'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-6321&redirect=T
model=changeRxns(model,'R02320','Phosphoenolpyruvate[c] + NDP[c] <=> Pyruvate[c] + Nucleoside triphosphate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02320'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.7.1.40&quickSearch=Quick+Search
model=setParam(model,'lb',{'R02322'},[0]);model=setParam(model,'ub',{'R02322'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NMNAMIDOHYDRO-RXN
model=setParam(model,'lb',{'R02323'},[0]);model=setParam(model,'ub',{'R02323'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-5841
model=setParam(model,'lb',{'R02325'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DCTP-DEAM-RXN
model=setParam(model,'lb',{'R02326'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DCDPKIN-RXN
model=setParam(model,'lb',{'R02328'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DTDPGLUCOSEPP-RXN
model=setParam(model,'lb',{'R02331'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DUDPKIN-RXN
%R02340 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-2381
model=setParam(model,'lb',{'R02401'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTARATE-SEMIALDEHYDE-DEHYDROGENASE-RXN
%R02404 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SUCCCOASYN-RXN
model=setParam(model,'lb',{'R02408'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=CYSTHIOCYS-RXN
model=changeRxns(model,'R02413','NADPH[c] + H+[c] + 3-Dehydroshikimate[c] <=> NADP+[c] + Shikimate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02413'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=SHIKIMATE-5-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R02427'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12246
model=changeRxns(model,'R02472','NADPH[c] + H+[c] + 2-Dehydropantoate[c] <=> NADP+[c] + (R)-Pantoate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02472'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2-DEHYDROPANTOATE-REDUCT-RXN&redirect=T
model=setParam(model,'lb',{'R02473'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PANTOATE-BETA-ALANINE-LIG-RXN&redirect=T
%R02484 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=URA-PHOSPH-RXN&redirect=T
model=setParam(model,'lb',{'R02487'},[0]);model=setParam(model,'ub',{'R02487'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTARYL-COA-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R02488'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTARYL-COA-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R02521'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=4-HYDROXYPHENYLPYRUVATE-DIOXYGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R02522'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=L-ARABINONATE-DEHYDRATASE-RXN&redirect=T
%R02527 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=D-LACTALDEHYDE-DEHYDROGENASE-RXN&redirect=T
%R02530 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLYOXI-RXN&redirect=T
model=setParam(model,'lb',{'R02539'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10819&redirect=T
%R02545 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.93&quickSearch=Quick+Search
model=setParam(model,'lb',{'R02549'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=AMINOBUTDEHYDROG-RXN&redirect=T
model=changeRxns(model,'R02550','Oxygen[c] + 2 H+[c] + 2 Reduced ferredoxin[c] + Toluene[c] <=> H2O[c] + 2 Oxidized ferredoxin[c] + Benzyl alcohol[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02550'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=TOLUENE-SIDE-CHAIN-MONOOXYGENASE-RXN
model=setParam(model,'lb',{'R02558'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PRUNASIN-BETA-GLUCOSIDASE-RXN&redirect=T
model=setParam(model,'lb',{'R02568'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8631&redirect=T
model=changeRxns(model,'R02569','CoA[c] + [Dihydrolipoyllysine-residue acetyltransferase] S-acetyldihydrolipoyllysine[c] <=> Acetyl-CoA[c] + Enzyme N6-(dihydrolipoyl)lysine[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02569'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-1133&redirect=T
model=changeRxns(model,'R02570','CoA[c] + [Dihydrolipoyllysine-residue succinyltransferase] S-succinyldihydrolipoyllysine[c] <=> Succinyl-CoA[c] + Enzyme N6-(dihydrolipoyl)lysine[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02570'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-1147&redirect=T
model=changeRxns(model,'R02571','CoA[c] + [Dihydrolipoyllysine-residue succinyltransferase] S-glutaryldihydrolipoyllysine[c] <=> Glutaryl-CoA[c] + Enzyme N6-(dihydrolipoyl)lysine[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02571'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.3.1.61&quickSearch=Quick+Search
                                    model=setParam(model,'lb',{'R02596'},[-1000]);model=setParam(model,'ub',{'R02596'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R02596'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R02601'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2-OXOPENT-4-ENOATE-HYDRATASE-RXN&redirect=T
model=setParam(model,'lb',{'R02602'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=4.1.1.77-RXN&redirect=T
model=changeRxns(model,'R02604','H2O[c] + 2-Hydroxymuconate semialdehyde[c] <=> Formate[c] + 2-Hydroxy-2,4-pentadienoate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02604'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3.7.1.9-RXN&redirect=T
%R02649 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ACETYLGLUTKIN-RXN&redirect=T
model=setParam(model,'lb',{'R02656'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GENTISATE-12-DIOXYGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R02668'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3-HYDROXY-KYNURENINASE-RXN&redirect=T
                                    model=setParam(model,'lb',{'R02670'},[-1000]);model=setParam(model,'ub',{'R02670'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R02670'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R02678'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10715&redirect=T
model=setParam(model,'lb',{'R02719'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=XMPXAN-RXN&redirect=T
model=setParam(model,'lb',{'R02720'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-1603&redirect=T
model=setParam(model,'lb',{'R02722'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TRYPSYN-RXN&redirect=T
model=setParam(model,'lb',{'R02734'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=SUCCDIAMINOPIMDESUCC-RXN&redirect=T
%R02735 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DIAMINOPIMEPIM-RXN&redirect=T
model=setParam(model,'lb',{'R02736'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLU6PDEHYDROG-RXN&redirect=T
%R02739 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLUCOSE-6-PHOSPHATE-1-EPIMERASE-RXN&redirect=T
%R02740 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PGLUCISOM-RXN
model=setParam(model,'lb',{'R02752'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLUCARDEHYDRA-RXN&redirect=T
model=setParam(model,'lb',{'R02762'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8527&redirect=T
model=changeRxns(model,'R02763','3-Carboxy-2-hydroxymuconate semialdehyde[c] <=> CO2[c] + 2-Hydroxymuconate semialdehyde[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02763'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15423
model=setParam(model,'lb',{'R02777'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DTDPDEHYRHAMREDUCT-RXN&redirect=T
model=setParam(model,'lb',{'R02783'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UDP-NACMURALA-GLU-LIG-RXN&redirect=T
model=setParam(model,'lb',{'R02788'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UDP-NACMURALGLDAPLIG-RXN&redirect=T
model=changeRxns(model,'R02804','2 H+[c] + 2 Ferrocytochrome c[c] + Nitrous oxide[c] <=> H2O[c] + 2 Ferricytochrome c[c] + Nitrogen[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02804'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12130
model=setParam(model,'lb',{'R02869'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=SPERMINE-SYNTHASE-RXN&redirect=T
model=changeRxns(model,'R02914','H2O[c] + Urocanate[c] <=> 4-Imidazolone-5-propanoate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                    
model=setParam(model,'lb',{'R02914'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UROCANATE-HYDRATASE-RXN&redirect=T
model=setParam(model,'lb',{'R02918'},[0]);model=setParam(model,'ub',{'R02918'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TYROSINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R02921'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.7.7.40-RXN&redirect=T
model=setParam(model,'lb',{'R02922'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=CREATININE-DEAMINASE-RXN&redirect=T
model=setParam(model,'lb',{'R02933'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8783&redirect=T
model=setParam(model,'lb',{'R02940'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14224&redirect=T
model=setParam(model,'lb',{'R02952'},[0]);model=setParam(model,'ub',{'R02952'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-802&redirect=T
model=setParam(model,'lb',{'R02957'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14225&redirect=T
model=setParam(model,'lb',{'R02964'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=NADH-DEHYDROGENASE-QUINONE-RXN&redirect=T
model=setParam(model,'lb',{'R02968'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=NAPHTHALENE-12-DIOXYGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R02971'},[0]);model=setParam(model,'ub',{'R02971'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PANTETHEINE-KINASE-RXN&redirect=T
                                    model=setParam(model,'lb',{'R02984'},[-1000]);model=setParam(model,'ub',{'R02984'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R02984'))),:)=[3];  % Reaction directionality unclear!
model=changeRxns(model,'R02985','H2O[c] + Amygdalin[c] <=> D-Glucose[c] + Prunasin[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R02985'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=AMYGDALIN-BETA-GLUCOSIDASE-RXN&redirect=T
model=setParam(model,'lb',{'R02986'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2-FUROATE--COA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R02987'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2-FUROYL-COA-DEHYDROGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R02990'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3-OXOADIPATE-COA-TRANSFERASE-RXN&redirect=T
model=setParam(model,'lb',{'R03005'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=NICONUCADENYLYLTRAN-RXN&redirect=T
%R03012 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=HISTOLDEHYD-RXN&redirect=T
model=setParam(model,'lb',{'R03018'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PANTOTHENATE-KIN-RXN&redirect=T
%R03026 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11667&redirect=T
model=setParam(model,'lb',{'R03033'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GALACTONDEHYDRAT-RXN&redirect=T
model=setParam(model,'lb',{'R03035'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PANTEPADENYLYLTRAN-RXN&redirect=T
model=setParam(model,'lb',{'R03038'},[0]);model=setParam(model,'ub',{'R03038'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ALANINE--TRNA-LIGASE-RXN&redirect=T
%R03045 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-6383&redirect=T
                                    model=setParam(model,'lb',{'R03050'},[-1000]);model=setParam(model,'ub',{'R03050'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03050'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R03051'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ACETOLACTREDUCTOISOM-RXN
                                    model=setParam(model,'lb',{'R03066'},[-1000]);model=setParam(model,'ub',{'R03066'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03066'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R03067'},[-1000]);model=setParam(model,'ub',{'R03067'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03067'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R03083'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3-DEHYDROQUINATE-SYNTHASE-RXN&redirect=T
%R03084 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3-DEHYDROQUINATE-DEHYDRATASE-RXN&redirect=T
model=setParam(model,'lb',{'R03105'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MERCAPYSTRANS-RXN
%R03132 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=SULFOCYS-RXN&redirect=T
model=setParam(model,'lb',{'R03140'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-961&redirect=T
model=setParam(model,'lb',{'R03165'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UROGENIIISYN-RXN&redirect=T
model=setParam(model,'lb',{'R03182'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DETHIOBIOTIN-SYN-RXN&redirect=T
model=changeRxns(model,'R03191','NADH[c] + H+[c] + UDP-N-acetyl-3-(1-carboxyvinyl)-D-glucosamine[c] <=> NAD+[c] + UDP-N-acetylmuramate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03191'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.98&quickSearch=Quick+Search
model=changeRxns(model,'R03192','NADPH[c] + H+[c] + UDP-N-acetyl-3-(1-carboxyvinyl)-D-glucosamine[c] <=> NADP+[c] + UDP-N-acetylmuramate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03192'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.98&quickSearch=Quick+Search
model=setParam(model,'lb',{'R03193'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UDP-NACMUR-ALA-LIG-RXN&redirect=T
model=setParam(model,'lb',{'R03194'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13403
model=setParam(model,'lb',{'R03197'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UROGENDECARBOX-RXN&redirect=T
model=setParam(model,'lb',{'R03210'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=7KAPSYN-RXN&redirect=T
model=setParam(model,'lb',{'R03220'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-1461&redirect=T
model=setParam(model,'lb',{'R03223'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=THI-P-SYN-RXN&redirect=T
%R03231 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DAPASYN-RXN&redirect=T
model=setParam(model,'lb',{'R03236'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TAGAKIN-RXN&redirect=T
model=setParam(model,'lb',{'R03237'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TAGAKIN-RXN&redirect=T
model=setParam(model,'lb',{'R03238'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TAGAKIN-RXN&redirect=T
model=setParam(model,'lb',{'R03239'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TAGAKIN-RXN&redirect=T
%R03243 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=HISTAMINOTRANS-RXN&redirect=T
model=changeRxns(model,'R03254','H2O[c] + Phosphoenolpyruvate[c] + D-Arabinose 5-phosphate[c] <=> Orthophosphate[c] + 3-Deoxy-D-manno-octulosonate 8-phosphate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03254'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=KDO-8PSYNTH-RXN&redirect=T
model=setParam(model,'lb',{'R03269'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=P-PANTOCYSDECARB-RXN&redirect=T
model=setParam(model,'lb',{'R03270'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12508&redirect=T
%R03276 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=5.1.2.3-RXN&redirect=T
model=setParam(model,'lb',{'R03283'},[0]);model=setParam(model,'ub',{'R03283'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.2.1.47-RXN&redirect=T
model=changeRxns(model,'R03291','NADH[c] + H+[c] + L-1-Pyrroline-3-hydroxy-5-carboxylate[c] <=> NAD+[c] + Hydroxyproline[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03291'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN66-546
model=changeRxns(model,'R03293','NADPH[c] + H+[c] + L-1-Pyrroline-3-hydroxy-5-carboxylate[c] <=> NADP+[c] + Hydroxyproline[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03293'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN66-546
model=setParam(model,'lb',{'R03303'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.13.11.15-RXN&redirect=T
model=changeRxns(model,'R03313','NADPH[c] + H+[c] + L-Glutamyl 5-phosphate[c] <=> NADP+[c] + Orthophosphate[c] + L-Glutamate 5-semialdehyde[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03313'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLUTSEMIALDEHYDROG-RXN&redirect=T
                                    model=setParam(model,'lb',{'R03316'},[-1000]);model=setParam(model,'ub',{'R03316'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03316'))),:)=[3];  % Reaction directionality unclear!
%R03321 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PGLUCISOM-RXN
model=setParam(model,'lb',{'R03346'},[0]);model=setParam(model,'ub',{'R03346'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14227&redirect=T
model=changeRxns(model,'R03348','5-Phospho-alpha-D-ribose 1-diphosphate[c] + Quinolinate[c] <=> CO2[c] + Diphosphate[c] + Nicotinate D-ribonucleotide[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03348'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=QUINOPRIBOTRANS-RXN&redirect=T
model=setParam(model,'lb',{'R03350'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=KDO-8PPHOSPHAT-RXN&redirect=T
model=setParam(model,'lb',{'R03351'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=CPM-KDOSYNTH-RXN&redirect=T
model=setParam(model,'lb',{'R03367'},[0]);model=setParam(model,'ub',{'R03367'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-21078
model=setParam(model,'lb',{'R03383'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12359
model=setParam(model,'lb',{'R03387'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DEHYDDEOXGALACTKIN-RXN&redirect=T 
model=changeRxns(model,'R03391','NADH[c] + H+[c] + CDP-4-dehydro-6-deoxy-D-glucose[c] <=> H2O[c] + NAD+[c] + CDP-4-dehydro-3,6-dideoxy-D-glucose[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03391'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.17.1.1-RXN
model=changeRxns(model,'R03392','NADPH[c] + H+[c] + CDP-4-dehydro-6-deoxy-D-glucose[c] <=> H2O[c] + NADP+[c] + CDP-4-dehydro-3,6-dideoxy-D-glucose[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03392'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.17.1.1-RXN&redirect=T
model=setParam(model,'lb',{'R03396'},[0]);model=setParam(model,'ub',{'R03396'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GDP-6-DEOXY-D-TALOSE-4-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R03398'},[0]);model=setParam(model,'ub',{'R03398'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GDP-6-DEOXY-D-TALOSE-4-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R03409'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PPPGPPHYDRO-RXN&redirect=T
%R03425 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GCVP-RXN&redirect=T
%R03443 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=N-ACETYLGLUTPREDUCT-RXN&redirect=T
model=setParam(model,'lb',{'R03457'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=IMIDPHOSDEHYD-RXN&redirect=T
%model=changeRxns(model,'R03458','NADPH[c] + H+[c] + 5-Amino-6-(5'-phosphoribosylamino)uracil[c] <=> NADP+[c] + 5-Amino-6-(5'-phospho-D-ribitylamino)uracil[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        

                                                                           % The changeRxns-command above did not work, therefore reactants are turned into products and
                                                                           % vice versa in the following way:
%constructEquations(model,(ismember(model.rxns,'R03458')))
Irxn=ismember(model.rxns,'R03458');                                         
Imet=ismember(model.mets,'C00005');        % NADPH[c]                                
model.S(Imet,Irxn)=[-1];                                                   % The coefficient for NADPH[c] is changed from 1 to -1.
Imet=ismember(model.mets,'C00080');        % H+[c]
model.S(Imet,Irxn)=[-1];                                                   % The coefficient for H+[c] is changed from 1 to -1.
Imet=ismember(model.mets,'C01268');        % 5-Amino-6-(5'-phosphoribosylamino)uracil[c]
model.S(Imet,Irxn)=[-1];                                                   % The coefficient for 5-Amino-6-(5'-phosphoribosylamino)uracil[c] is changed from 1 to -1.
Imet=ismember(model.mets,'C00006');        % NADP+[c]
model.S(Imet,Irxn)=[1];                                                    % The coefficient for NADP+[c] is changed from -1 to 1.
Imet=ismember(model.mets,'C04454');        % 5-Amino-6-(5'-phospho-D-ribitylamino)uracil[c]
model.S(Imet,Irxn)=[1];                                                    % The coefficient for 5-Amino-6-(5'-phospho-D-ribitylamino)uracil[c] is changed from -1 to 1.
%constructEquations(model,(ismember(model.rxns,'R03458')))

model=setParam(model,'lb',{'R03458'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RIBOFLAVINSYNREDUC-RXN&redirect=T
model=setParam(model,'lb',{'R03459'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RIBOFLAVINSYNDEAM-RXN&redirect=T
%R03460 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.5.1.19-RXN&redirect=T
model=setParam(model,'lb',{'R03470'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=4-CARBOXYMUCONOLACTONE-DECARBOXYLASE-RXN&redirect=T
model=setParam(model,'lb',{'R03471'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=OHMETPYRKIN-RXN&redirect=T
model=setParam(model,'lb',{'R03472'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PYRIMSYN1-RXN
model=setParam(model,'lb',{'R03493'},[0]);model=setParam(model,'ub',{'R03493'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ISOHEXENYLGLUTACONYL-COA-HYDRATASE-RXN&redirect=T
model=setParam(model,'lb',{'R03494'},[0]);model=setParam(model,'ub',{'R03494'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GERANOYL-COA-CARBOXYLASE-RXN&redirect=T
model=setParam(model,'lb',{'R03503'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=H2PTERIDINEPYROPHOSPHOKIN-RXN&redirect=T
model=setParam(model,'lb',{'R03504'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=H2NEOPTERINALDOL-RXN&redirect=T
model=setParam(model,'lb',{'R03508'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=IGPSYN-RXN&redirect=T
model=setParam(model,'lb',{'R03509'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PRAISOM-RXN&redirect=T
model=setParam(model,'lb',{'R03522'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GSHTRAN-RXN&redirect=T
model=setParam(model,'lb',{'R03527'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.2.1.21&quickSearch=Quick+Search
%R03530 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14228&redirect=T
model=setParam(model,'lb',{'R03531'},[0]);model=setParam(model,'ub',{'R03531'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-1602&redirect=T
model=setParam(model,'lb',{'R03549'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11861&redirect=T
model=setParam(model,'lb',{'R03550'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GALLATE-DIOXYGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R03560'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=TOLUENE-2-MONOOXYGENASE-RXN
model=setParam(model,'lb',{'R03562'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.13.-&quickSearch=Quick+Search
model=setParam(model,'lb',{'R03566'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3-5-XYLENOL-METHYLHYDROXYLASE-RXN
model=setParam(model,'lb',{'R03595'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.7.9.3-RXN&redirect=T
                                    model=setParam(model,'lb',{'R03596'},[-1000]);model=setParam(model,'ub',{'R03596'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03596'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R03599'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=SELENOCYSTEINE-LYASE-RXN&redirect=T
model=setParam(model,'lb',{'R03601'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12726
model=setParam(model,'lb',{'R03608'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=O-CRESOL-METHYLCATECHOL-RXN
model=setParam(model,'lb',{'R03643'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-19671
model=setParam(model,'lb',{'R03646'},[0]);model=setParam(model,'ub',{'R03646'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ARGININE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03650'},[0]);model=setParam(model,'ub',{'R03650'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=CYSTEINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03652'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLUTAMINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03654'},[0]);model=setParam(model,'ub',{'R03654'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLYCINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03655'},[0]);model=setParam(model,'ub',{'R03655'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=HISTIDINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03656'},[0]);model=setParam(model,'ub',{'R03656'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ISOLEUCINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03657'},[0]);model=setParam(model,'ub',{'R03657'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LEUCINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03658'},[0]);model=setParam(model,'ub',{'R03658'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LYSINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03659'},[0]);model=setParam(model,'ub',{'R03659'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=METHIONINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03660'},[0]);model=setParam(model,'ub',{'R03660'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PHENYLALANINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03661'},[0]);model=setParam(model,'ub',{'R03661'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PROLINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03662'},[0]);model=setParam(model,'ub',{'R03662'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=SERINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03663'},[0]);model=setParam(model,'ub',{'R03663'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=THREONINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03664'},[0]);model=setParam(model,'ub',{'R03664'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TRYPTOPHAN--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03665'},[0]);model=setParam(model,'ub',{'R03665'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=VALINE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R03751'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=6-HEXANOLIDE-HYDROLYSIS-RXN&redirect=T
model=changeRxns(model,'R03778','CoA[c] + 3-Oxodecanoyl-CoA[c] <=> Acetyl-CoA[c] + Octanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03778'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-13617&redirect=T
model=setParam(model,'lb',{'R03789'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=QUEUOSINE-TRNA-RIBOSYLTRANSFERASE-RXN&redirect=T
model=changeRxns(model,'R03791','(R)-Mandelate[c] <=> (S)-Mandelate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03791'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TAGAKIN-RXN&redirect=T
model=setParam(model,'lb',{'R03813'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LEUCYLTRANSFERASE-RXN&redirect=T
model=setParam(model,'lb',{'R03815'},[-1000]);model=setParam(model,'ub',{'R03815'},[1000]); %https://biocyc.org/META/substring-search?type=NIL&object=1.8.1.4&quickSearch=Quick+Search                              Irreversible in the original FDR.
model=setParam(model,'lb',{'R03816'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-19671           
model=changeRxns(model,'R03858','CoA[c] + 3-Oxotetradecanoyl-CoA[c] <=> Acetyl-CoA[c] + Lauroyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03858'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14268&redirect=T
model=setParam(model,'lb',{'R03867'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN66-336
                                    model=setParam(model,'lb',{'R03869'},[-1000]);model=setParam(model,'ub',{'R03869'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03869'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R03889'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.2.1.32-RXN&redirect=T
model=setParam(model,'lb',{'R03893'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=CARBOXYMETHYLENEBUTENOLIDASE-RXN&redirect=T
%R03896 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=R-2-METHYLMALATE-DEHYDRATASE-RXN&redirect=T
model=setParam(model,'lb',{'R03898'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-7744&redirect=T
model=changeRxns(model,'R03905','H2O[c] + ATP[c] + L-Glutamine[c] + L-Glutamyl-tRNA(Gln)[c] <=> ADP[c] + Orthophosphate[c] + L-Glutamate[c] + Glutaminyl-tRNA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                        
model=setParam(model,'lb',{'R03905'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=6.3.5.7-RXN&redirect=T
                                    model=setParam(model,'lb',{'R03916'},[-1000]);model=setParam(model,'ub',{'R03916'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03916'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R03919'},[-1000]);model=setParam(model,'ub',{'R03919'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03919'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R03920'},[0]);model=setParam(model,'ub',{'R03920'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUCOKIN-RXN
                                    model=setParam(model,'lb',{'R03936'},[-1000]);model=setParam(model,'ub',{'R03936'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03936'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R03940'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=METHIONYL-TRNA-FORMYLTRANSFERASE-RXN&redirect=T
model=setParam(model,'lb',{'R03948'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.1.1.130-RXN&redirect=T
%R03966 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8529&redirect=T
%R03968 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3-ISOPROPYLMALISOM-RXN&redirect=T
                                    model=setParam(model,'lb',{'R03970'},[-1000]);model=setParam(model,'ub',{'R03970'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03970'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R03971'},[-1000]);model=setParam(model,'ub',{'R03971'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03971'))),:)=[3];  % Reaction directionality unclear!                                    
model=changeRxns(model,'R03991','CoA[c] + 3-Oxopalmitoyl-CoA[c] <=> Acetyl-CoA[c] + Tetradecanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                            
model=setParam(model,'lb',{'R03991'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.3.1.155-RXN&redirect=T
model=setParam(model,'lb',{'R03998'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.13.40&quickSearch=Quick+Search
model=setParam(model,'lb',{'R03999'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.13.40&quickSearch=Quick+Search
%R04001 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8991
                                    model=setParam(model,'lb',{'R04007'},[-1000]);model=setParam(model,'ub',{'R04007'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04007'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R04035'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=HISTPRATPHYD-RXN&redirect=T
model=setParam(model,'lb',{'R04037'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=HISTCYCLOHYD-RXN&redirect=T
model=setParam(model,'lb',{'R04065'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10089&redirect=T
model=setParam(model,'lb',{'R04089'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=CATECHOL-2-3-DIOXYGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R04095'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.8.7&quickSearch=Quick+Search
model=setParam(model,'lb',{'R04109'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLUTRNAREDUCT-RXN&redirect=T
model=setParam(model,'lb',{'R04112'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14267&redirect=T
%R04125 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GCVT-RXN
                                    model=setParam(model,'lb',{'R04132'},[-1000]);model=setParam(model,'ub',{'R04132'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04132'))),:)=[3];  % Reaction directionality unclear!
model=changeRxns(model,'R04134','2-Hydroxyhepta-2,4-dienedioate[c] <=> 2-Oxohept-3-enedioate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04134'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN1K-87&redirect=T
model=changeRxns(model,'R04137','H2O[c] + 3-Methylcrotonyl-CoA[c] <=> 3-Hydroxyisovaleryl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04137'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14266&redirect=T
model=setParam(model,'lb',{'R04138'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=METHYLCROTONYL-COA-CARBOXYLASE-RXN&redirect=T
model=setParam(model,'lb',{'R04144'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLYRIBONUCSYN-RXN&redirect=T
model=setParam(model,'lb',{'R04148'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DMBPPRIBOSYLTRANS-RXN&redirect=T
model=changeRxns(model,'R04161','(R)-4-Hydroxymandelate[c] <=> (S)-4-Hydroxymandelate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04161'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=HYDROXY-MANDELATE-RACEMASE-RXN&redirect=T
model=changeRxns(model,'R04170','H2O[c] + 2-trans-Dodecenoyl-CoA[c] <=> (S)-3-Hydroxydodecanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04170'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17797
%R04173 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PSERTRANSAM-RXN&redirect=T
%R04187 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14265&redirect=T
%R04188 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.6.1.22-RXN&redirect=T
model=changeRxns(model,'R04198','NADH[c] + H+[c] + (2S,4S)-4-Hydroxy-2,3,4,5-tetrahydrodipicolinate[c] <=> H2O[c] + NAD+[c] + 2,3,4,5-Tetrahydrodipicolinate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04198'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14014&redirect=T
model=changeRxns(model,'R04199','NADPH[c] + H+[c] + (2S,4S)-4-Hydroxy-2,3,4,5-tetrahydrodipicolinate[c] <=> H2O[c] + NADP+[c] + 2,3,4,5-Tetrahydrodipicolinate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04199'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14014&redirect=T
%R04203 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.1.1.178-RXN&redirect=T
%R04204 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TIGLYLCOA-HYDROXY-RXN&redirect=T
model=setParam(model,'lb',{'R04208'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=AIRS-RXN&redirect=T
model=changeRxns(model,'R04212','H2O[c] + ATP[c] + L-Glutamine[c] + L-Aspartyl-tRNA(Asn)[c] <=> ADP[c] + Orthophosphate[c] + L-Glutamate[c] + L-Asparaginyl-tRNA(Asn)[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04212'},[0]);model=setParam(model,'ub',{'R04212'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=6.3.5.6-RXN&redirect=T
model=setParam(model,'lb',{'R04218'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXNARA-8002
%R04224 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=METHYLACYLYLCOA-HYDROXY-RXN&redirect=T
model=setParam(model,'lb',{'R04231'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=P-PANTOCYSLIG-RXN&redirect=T
model=setParam(model,'lb',{'R04277'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-2462&redirect=T
model=setParam(model,'lb',{'R04278'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.2.1.45-RXN&redirect=T
model=setParam(model,'lb',{'R04279'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.2.1.45-RXN&redirect=T
model=setParam(model,'lb',{'R04280'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=METHYLGALLATE-RXN&redirect=T
model=changeRxns(model,'R04292','Glycerone phosphate[c] + Iminoaspartate[c] <=> 2 H2O[c] + Orthophosphate[c] + Quinolinate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04292'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=QUINOLINATE-SYNTHA-RXN&redirect=T
%R04325 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GART-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04326'},[-1000]);model=setParam(model,'ub',{'R04326'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04326'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R04355'},[-1000]);model=setParam(model,'ub',{'R04355'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04355'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R04365'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TETHYDPICSUCC-RXN&redirect=T
%model=changeRxns(model,'R04378','5-Phospho-alpha-D-ribose 1-diphosphate[c] + 5-Amino-4-imidazolecarboxyamide[c] <=> Diphosphate[c] + 1-(5'-Phosphoribosyl)-5-amino-4-imidazolecarboxamide[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                

                                                                           % The changeRxns-command above did not work, therefore reactants are turned into products and
                                                                           % vice versa in the following way:
%constructEquations(model,(ismember(model.rxns,'R04378')))
Irxn=ismember(model.rxns,'R04378');                                         
Imet=ismember(model.mets,'C00119');        % 5-Phospho-alpha-D-ribose 1-diphosphate[c]                               
model.S(Imet,Irxn)=[-1];                                                   % The coefficient for NADPH[c] is changed from 1 to -1.
Imet=ismember(model.mets,'C04051');        % 5-Amino-4-imidazolecarboxyamide[c]
model.S(Imet,Irxn)=[-1];                                                   % The coefficient for H+[c] is changed from 1 to -1.
Imet=ismember(model.mets,'C00013');        % Diphosphate[c]
model.S(Imet,Irxn)=[1];                                                   % The coefficient for Diphosphate[c] is changed from -1 to 1.
Imet=ismember(model.mets,'C04677');        % 1-(5'-Phosphoribosyl)-5-amino-4-imidazolecarboxamide[c]
model.S(Imet,Irxn)=[1];                                                    % The coefficient for 1-(5'-Phosphoribosyl)-5-amino-4-imidazolecarboxamide[c] is changed from -1 to 1.
%constructEquations(model,(ismember(model.rxns,'R04378')))

model=setParam(model,'lb',{'R04378'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14270&redirect=T
model=setParam(model,'lb',{'R04379'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=5.3.3.10-RXN&redirect=T
model=setParam(model,'lb',{'R04380'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=4.1.1.68&quickSearch=Quick+Search
model=setParam(model,'lb',{'R04385'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=BIOTIN-CARBOXYL-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04386'},[-1000]);model=setParam(model,'ub',{'R04386'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04386'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R04391'},[-1000]);model=setParam(model,'ub',{'R04391'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04391'))),:)=[3];  % Reaction directionality unclear!
%R04405 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=HOMOCYSMET-RXN&redirect=T
model=setParam(model,'lb',{'R04418'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=CHMS-DEHYDROGENASE-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04419'},[-1000]);model=setParam(model,'ub',{'R04419'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04419'))),:)=[3];  % Reaction directionality unclear!
%R04425 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=4.2.1.99-RXN&redirect=T
model=setParam(model,'lb',{'R04426'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3-ISOPROPYLMALDEHYDROG-RXN&redirect=T
model=setParam(model,'lb',{'R04428'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=4.2.1.58-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04429'},[-1000]);model=setParam(model,'ub',{'R04429'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04429'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R04430'},[-1000]);model=setParam(model,'ub',{'R04430'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04430'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R04439'},[0]);model=setParam(model,'ub',{'R04439'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=R04439&quickSearch=Quick+Search
model=changeRxns(model,'R04440','NADPH[c] + H+[c] + 3-Hydroxy-3-methyl-2-oxobutanoic acid[c] <=> NADP+[c] + (R)-2,3-Dihydroxy-3-methylbutanoate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04440'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.86&quickSearch=Quick+Search
model=setParam(model,'lb',{'R04441'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DIHYDROXYISOVALDEHYDRAT-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04444'},[-1000]);model=setParam(model,'ub',{'R04444'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04444'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R04445'},[-1000]);model=setParam(model,'ub',{'R04445'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04445'))),:)=[3];  % Reaction directionality unclear!                                     
model=setParam(model,'lb',{'R04457'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=LUMAZINESYN-RXN
model=setParam(model,'lb',{'R04463'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=FGAMSYN-RXN&redirect=T
model=changeRxns(model,'R04475','L-Glutamate[c] + N-Succinyl-2-L-amino-6-oxoheptanedioate[c] <=> 2-Oxoglutarate[c] + N-Succinyl-LL-2,6-diaminoheptanedioate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04475'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=SUCCINYLDIAMINOPIMTRANS-RXN&redirect=T
%R04478 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-2463&redirect=T
                                    model=setParam(model,'lb',{'R04482'},[-1000]);model=setParam(model,'ub',{'R04482'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04482'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R04506'},[-1000]);model=setParam(model,'ub',{'R04506'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04506'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R04509'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PYRIMSYN3-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04533'},[-1000]);model=setParam(model,'ub',{'R04533'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04533'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R04534'},[-1000]);model=setParam(model,'ub',{'R04534'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04534'))),:)=[3];  % Reaction directionality unclear!                                    
                                    model=setParam(model,'lb',{'R04535'},[-1000]);model=setParam(model,'ub',{'R04535'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04535'))),:)=[3];  % Reaction directionality unclear!                                    
                                    model=setParam(model,'lb',{'R04536'},[-1000]);model=setParam(model,'ub',{'R04536'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04536'))),:)=[3];  % Reaction directionality unclear!                                    
                                    model=setParam(model,'lb',{'R04537'},[-1000]);model=setParam(model,'ub',{'R04537'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04537'))),:)=[3];  % Reaction directionality unclear!                                    
                                    model=setParam(model,'lb',{'R04543'},[-1000]);model=setParam(model,'ub',{'R04543'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04543'))),:)=[3];  % Reaction directionality unclear!                                    
                                    model=setParam(model,'lb',{'R04544'},[-1000]);model=setParam(model,'ub',{'R04544'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04544'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R04549'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LIPIDXSYNTHESIS-RXN&redirect=T
model=setParam(model,'lb',{'R04550'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UDPHYDROXYMYRGLUCOSAMNACETYLTRANS-RXN&redirect=T
model=setParam(model,'lb',{'R04558'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLUTAMIDOTRANS-RXN&redirect=T
%R04559 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=AICARSYN-RXN&redirect=T
%R04560 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=AICARTRANSFORM-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04566'},[-1000]);model=setParam(model,'ub',{'R04566'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04566'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R04567'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UDPNACETYLGLUCOSAMACYLTRANS-RXN&redirect=T
model=setParam(model,'lb',{'R04568'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN66-637
model=setParam(model,'lb',{'R04573'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=6.3.2.10-RXN&redirect=T
model=setParam(model,'lb',{'R04587'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UDPACYLGLCNACDEACETYL-RXN&redirect=T
model=setParam(model,'lb',{'R04591'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=SAICARSYN-RXN&redirect=T
model=setParam(model,'lb',{'R04594'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16788
model=setParam(model,'lb',{'R04606'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LIPIDADISACCHARIDESYNTH-RXN&redirect=T
model=setParam(model,'lb',{'R04617'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UDP-NACMURALGLDAPAALIG-RXN&redirect=T
model=setParam(model,'lb',{'R04638'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=H2NEOPTERINP3PYROPHOSPHOHYDRO-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04640'},[-1000]);model=setParam(model,'ub',{'R04640'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04640'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R04657'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TETRAACYLDISACC4KIN-RXN&redirect=T
model=setParam(model,'lb',{'R04658'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=KDOTRANS-RXN&redirect=T
model=setParam(model,'lb',{'R04666'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.5.1.6&quickSearch=Quick+Search
                                    model=setParam(model,'lb',{'R04672'},[-1000]);model=setParam(model,'ub',{'R04672'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04672'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R04673'},[-1000]);model=setParam(model,'ub',{'R04673'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04673'))),:)=[3];  % Reaction directionality unclear!                                     
model=setParam(model,'lb',{'R04710'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=TDCEACT-RXN
model=changeRxns(model,'R04724','NADH[c] + H+[c] + trans-Dodec-2-enoyl-[acp][c] <=> NAD+[c] + Dodecanoyl-[acyl-carrier protein][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04724'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.9&quickSearch=Quick+Search
model=changeRxns(model,'R04725','NADPH[c] + H+[c] + trans-Dodec-2-enoyl-[acp][c] <=> NADP+[c] + Dodecanoyl-[acyl-carrier protein][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04725'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.10&quickSearch=Quick+Search
model=setParam(model,'lb',{'R04726'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9535
model=setParam(model,'lb',{'R04734'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14206
%R04737 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14271&redirect=T
model=changeRxns(model,'R04738','H2O[c] + trans-Hexadec-2-enoyl-CoA[c] <=> (S)-3-Hydroxyhexadecanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04738'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14272&redirect=T
model=setParam(model,'lb',{'R04739'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12507&redirect=T
model=changeRxns(model,'R04740','H2O[c] + trans-Tetradec-2-enoyl-CoA[c] <=> (S)-3-Hydroxytetradecanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                
model=setParam(model,'lb',{'R04740'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14273&redirect=T
                                    model=setParam(model,'lb',{'R04741'},[-1000]);model=setParam(model,'ub',{'R04741'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04741'))),:)=[3];  % Reaction directionality unclear!
model=changeRxns(model,'R04742','CoA[c] + 3-Oxododecanoyl-CoA[c] <=> Acetyl-CoA[c] + Decanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                   
model=setParam(model,'lb',{'R04742'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14274&redirect=T
model=setParam(model,'lb',{'R04743'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12490&redirect=T
model=changeRxns(model,'R04744','H2O[c] + trans-Dec-2-enoyl-CoA[c] <=> (S)-Hydroxydecanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                   
model=setParam(model,'lb',{'R04744'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-13616&redirect=T
model=setParam(model,'lb',{'R04745'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14275&redirect=T
model=changeRxns(model,'R04746','H2O[c] + trans-Oct-2-enoyl-CoA[c] <=> (S)-3-Hydroxyoctanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                   
model=setParam(model,'lb',{'R04746'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-20678
model=changeRxns(model,'R04747','CoA[c] + 3-Oxooctanoyl-CoA[c] <=> Acetyl-CoA[c] + Hexanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                   
model=setParam(model,'lb',{'R04747'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14277&redirect=T
%R04748 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12570&redirect=T
model=changeRxns(model,'R04749','H2O[c] + trans-Hex-2-enoyl-CoA[c] <=> (S)-Hydroxyhexanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                   
model=setParam(model,'lb',{'R04749'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12567&redirect=T
model=setParam(model,'lb',{'R04771'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=S-ADENMETSYN-RXN
model=setParam(model,'lb',{'R04773'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=6.1.1.10&quickSearch=Quick+Search
model=setParam(model,'lb',{'R04779'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=6PFRUCTPHOS-RXN
model=setParam(model,'lb',{'R04780'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=F16BDEPHOS-RXN
model=setParam(model,'lb',{'R04858'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DNA-CYTOSINE-5--METHYLTRANSFERASE-RXN
                                    model=setParam(model,'lb',{'R04859'},[-1000]);model=setParam(model,'ub',{'R04859'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04859'))),:)=[3];  % Reaction directionality unclear!
model=changeRxns(model,'R04880','NADH[c] + H+[c] + 3,4-Dihydroxymandelaldehyde[c] <=> NAD+[c] + 3,4-Dihydroxyphenylethyleneglycol[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04880'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10911&redirect=T
model=setParam(model,'lb',{'R04903'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10780&redirect=T
model=setParam(model,'lb',{'R04911'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17150
model=setParam(model,'lb',{'R04929'},[0]);model=setParam(model,'ub',{'R04929'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12720
model=setParam(model,'lb',{'R04935'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GAMMA-GLUTAMYLTRANSFERASE-RXN
model=setParam(model,'lb',{'R04936'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ADENOSYLHOMOCYSTEINASE-RXN
model=setParam(model,'lb',{'R04941'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12729
model=setParam(model,'lb',{'R04949'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.2.1.21&quickSearch=Quick+Search
model=setParam(model,'lb',{'R04951'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.4.13.-&quickSearch=Quick+Search
model=setParam(model,'lb',{'R04952'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.3.1.41&quickSearch=Quick+Search
model=changeRxns(model,'R04953','NADPH[c] + H+[c] + 3-Oxohexanoyl-[acp][c] <=> NADP+[c] + (R)-3-Hydroxyhexanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04953'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9518
model=setParam(model,'lb',{'R04954'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN66-620
model=changeRxns(model,'R04955','NADH[c] + H+[c] + trans-Hex-2-enoyl-[acp][c] <=> NAD+[c] + Hexanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04955'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9658
model=changeRxns(model,'R04956','NADPH[c] + H+[c] + trans-Hex-2-enoyl-[acp][c] <=> NADP+[c] + Hexanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04956'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9521
model=setParam(model,'lb',{'R04957'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9523
model=changeRxns(model,'R04958','NADH[c] + H+[c] + trans-Oct-2-enoyl-[acp][c] <=> NAD+[c] + Octanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04958'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9659
model=changeRxns(model,'R04959','NADPH[c] + H+[c] + trans-Oct-2-enoyl-[acp][c] <=> NADP+[c] + Octanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04959'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9526
model=setParam(model,'lb',{'R04960'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9527
model=changeRxns(model,'R04961','NADH[c] + H+[c] + trans-Dec-2-enoyl-[acp][c] <=> NAD+[c] + Decanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04961'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9660
model=changeRxns(model,'R04962','NADPH[c] + H+[c] + trans-Dec-2-enoyl-[acp][c] <=> NADP+[c] + Decanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04962'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9530
model=setParam(model,'lb',{'R04963'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9531
model=changeRxns(model,'R04964','NADPH[c] + H+[c] + 3-Oxododecanoyl-[acp][c] <=> NADP+[c] + (R)-3-Hydroxydodecanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04964'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9532
model=setParam(model,'lb',{'R04965'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN66-632
model=changeRxns(model,'R04966','NADH[c] + H+[c] + trans-Tetradec-2-enoyl-[acp][c] <=> NAD+[c] + Tetradecanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04966'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9662
model=changeRxns(model,'R04967','NADPH[c] + H+[c] + trans-Tetradec-2-enoyl-[acp][c] <=> NADP+[c] + Tetradecanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04967'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=ENOYL-ACP-REDUCT-NADPH-RXN
model=setParam(model,'lb',{'R04968'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.3.1.41&quickSearch=Quick+Search
model=changeRxns(model,'R04969','NADH[c] + H+[c] + trans-Hexadec-2-enoyl-[acp][c] <=> NAD+[c] + Hexadecanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04969'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.9&quickSearch=Quick+Search
model=changeRxns(model,'R04970','NADPH[c] + H+[c] + trans-Hexadec-2-enoyl-[acp][c] <=> NADP+[c] + Hexadecanoyl-[acp][c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R04970'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.10&quickSearch=Quick+Search
model=setParam(model,'lb',{'R04972'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10642&redirect=T
model=setParam(model,'lb',{'R04984'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.99.60&quickSearch=Quick+Search
model=setParam(model,'lb',{'R04986'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3-OCTAPRENYL-4-OHBENZOATE-DECARBOX-RXN&redirect=T
model=setParam(model,'lb',{'R04988'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2-OCTAPRENYL-6-OHPHENOL-METHY-RXN&redirect=T
model=setParam(model,'lb',{'R04990'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN3O-54&redirect=T 
model=setParam(model,'lb',{'R04998'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8036
model=setParam(model,'lb',{'R05032'},[-1000]);model=setParam(model,'ub',{'R05032'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=NACGLCTRANS-RXN&redirect=T                                   Irreversible in the original FDR.
                                    model=setParam(model,'lb',{'R05046'},[-1000]);model=setParam(model,'ub',{'R05046'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05046'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R05050'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-37&redirect=T
                                    model=setParam(model,'lb',{'R05051'},[-1000]);model=setParam(model,'ub',{'R05051'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05051'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R05066'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3-HYDROXYISOBUTYRATE-DEHYDROGENASE-RXN
model=changeRxns(model,'R05068','NADPH[c] + H+[c] + (R)-3-Hydroxy-3-methyl-2-oxopentanoate[c] <=> NADP+[c] + (R)-2,3-Dihydroxy-3-methylpentanoate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R05068'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.86&quickSearch=Quick+Search
%R05069 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14106
model=setParam(model,'lb',{'R05070'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DIHYDROXYMETVALDEHYDRAT-RXN&redirect=T
%R05071 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2-ACETOLACTATE-MUTASE-RXN&redirect=T
model=setParam(model,'lb',{'R05074'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=KDOTRANS2-RXN&redirect=T
%R05085 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PSERTRANSAMPYR-RXN&redirect=T
model=setParam(model,'lb',{'R05086'},[0]);model=setParam(model,'ub',{'R05086'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14125&redirect=T
model=setParam(model,'lb',{'R05118'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.1.1.75&quickSearch=Quick+Search
model=setParam(model,'lb',{'R05146'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LAUROYLACYLTRAN-RXN&redirect=T
model=setParam(model,'lb',{'R05149'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.1.1.132-RXN&redirect=T
%R05177 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=5.4.1.2-RXN&redirect=T
model=setParam(model,'lb',{'R05180'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.1.1.131-RXN&redirect=T
model=setParam(model,'lb',{'R05181'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.1.1.133-RXN&redirect=T
model=setParam(model,'lb',{'R05220'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=R344-RXN&redirect=T
model=setParam(model,'lb',{'R05221'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=COBINAMIDEKIN-RXN&redirect=T 
model=setParam(model,'lb',{'R05222'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=COBINPGUANYLYLTRANS-RXN&redirect=T
model=changeRxns(model,'R05223','alpha-Ribazole[c] + Adenosine-GDP-cobinamide[c] <=> GMP[c] + Cobamide coenzyme[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R05223'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=COBALAMINSYN-RXN&redirect=T
model=setParam(model,'lb',{'R05224'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=R341-RXN&redirect=T
model=setParam(model,'lb',{'R05225'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=R345-RXN&redirect=T
%R05233 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.1&quickSearch=Quick+Search
%R05234 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.1&quickSearch=Quick+Search
                                    model=setParam(model,'lb',{'R05237'},[-1000]);model=setParam(model,'ub',{'R05237'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05237'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R05238'},[-1000]);model=setParam(model,'ub',{'R05238'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05238'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R05248'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-666&redirect=T
model=setParam(model,'lb',{'R05274'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-10891
model=setParam(model,'lb',{'R05284'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DHLAXANAU-RXN&redirect=T
%R05286 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ALDXANAU-RXN&redirect=T
model=setParam(model,'lb',{'R05287'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DHLBXANAU-RXN&redirect=T
%R05292 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.3.1.67-RXN
%R05293 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.3.1.67-RXN
model=setParam(model,'lb',{'R05295'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11219&redirect=T
%R05298 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12735&redirect=T
%R05305 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.1.1.259-RXN&redirect=T
%R05309 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.3.1.68-RXN
%R05314 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.3.1.68-RXN&redirect=T
model=setParam(model,'lb',{'R05320'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-299&redirect=T
model=setParam(model,'lb',{'R05332'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.3.1.157-RXN&redirect=T
                                    model=setParam(model,'lb',{'R05353'},[-1000]);model=setParam(model,'ub',{'R05353'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05353'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R05362'},[-1000]);model=setParam(model,'ub',{'R05362'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05362'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R05367'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.8.1.5&quickSearch=Quick+Search
model=setParam(model,'lb',{'R05368'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.8.1.5&quickSearch=Quick+Search
model=setParam(model,'lb',{'R05369'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LINB1-RXN&redirect=T
model=setParam(model,'lb',{'R05370'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LINB2-RXN&redirect=T
model=setParam(model,'lb',{'R05374'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8841
model=setParam(model,'lb',{'R05377'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-667&redirect=T
%R05389 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=5.3.2.6&quickSearch=Quick+Search
model=setParam(model,'lb',{'R05404'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.13.11.2&quickSearch=Quick+Search
model=setParam(model,'lb',{'R05406'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.13.11.2&quickSearch=Quick+Search
model=setParam(model,'lb',{'R05422'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=FLUORENE-OXYGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R05423'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-2482&redirect=T
                                    model=setParam(model,'lb',{'R05424'},[-1000]);model=setParam(model,'ub',{'R05424'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05424'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R05425'},[-1000]);model=setParam(model,'ub',{'R05425'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05425'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R05426'},[-1000]);model=setParam(model,'ub',{'R05426'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05426'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R05427'},[-1000]);model=setParam(model,'ub',{'R05427'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05427'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R05506'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-2006&redirect=T
model=setParam(model,'lb',{'R05510'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.1.1.45&quickSearch=Quick+Search
model=setParam(model,'lb',{'R05511'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9857&redirect=T
model=setParam(model,'lb',{'R05553'},[0]);model=setParam(model,'ub',{'R05553'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ADCLY-RXN&redirect=T
model=changeRxns(model,'R05576','NADPH[c] + H+[c] + Acetoacetyl-CoA[c] <=> NADP+[c] + 3-Hydroxybutanoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R05576'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3-HYDROXYBUTYRYL-COA-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R05577'},[0]);model=setParam(model,'ub',{'R05577'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ASPARTATE--TRNA-LIGASE-RXN&redirect=T
model=setParam(model,'lb',{'R05578'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLURS-RXN&redirect=T
model=setParam(model,'lb',{'R05582'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=R267-RXN&redirect=T
%R05586 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8032&redirect=T
model=setParam(model,'lb',{'R05592'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=R268-RXN&redirect=T
                                    model=setParam(model,'lb',{'R05595'},[-1000]);model=setParam(model,'ub',{'R05595'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05595'))),:)=[3];  % Reaction directionality unclear!  
model=setParam(model,'lb',{'R05604'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=D-ARABINITOL-4-DEHYDROGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R05605'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=KDPGALDOL-RXN&redirect=T
model=setParam(model,'lb',{'R05608'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GALACTARDEHYDRA-RXN&redirect=T
model=setParam(model,'lb',{'R05614'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.1.1.64&quickSearch=Quick+Search
model=setParam(model,'lb',{'R05615'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=4OHBENZOATE-OCTAPRENYLTRANSFER-RXN&redirect=T
%R05619 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=R282-RXN&redirect=T
model=setParam(model,'lb',{'R05627'},[0]);model=setParam(model,'ub',{'R05627'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=UNDECAPRENYL-DIPHOSPHATASE-RXN&redirect=T
%R05629 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8975&redirect=T
model=setParam(model,'lb',{'R05630'},[0]);model=setParam(model,'ub',{'R05630'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PHOSNACMURPENTATRANS-RXN&redirect=T
model=setParam(model,'lb',{'R05632'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12763
model=setParam(model,'lb',{'R05633'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.7.7.60-RXN&redirect=T
model=setParam(model,'lb',{'R05634'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.7.1.148-RXN&redirect=T
model=setParam(model,'lb',{'R05636'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DXS-RXN&redirect=T
model=setParam(model,'lb',{'R05637'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-302&redirect=T
model=setParam(model,'lb',{'R05645'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-4301&redirect=T
model=setParam(model,'lb',{'R05647'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-4361&redirect=T
model=setParam(model,'lb',{'R05662'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8976&redirect=T
model=setParam(model,'lb',{'R05666'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=TOLUENE-4-MONOOXYGENASE-RXN
                                    model=setParam(model,'lb',{'R05681'},[-1000]);model=setParam(model,'ub',{'R05681'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05681'))),:)=[3];  % Reaction directionality unclear!
%R05688 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DXPREDISOM-RXN&redirect=T
model=changeRxns(model,'R05706','NADPH[c] + FMN[c] + H+[c] <=> NADP+[c] + Reduced FMN[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R05706'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=FMNREDUCT-RXN&redirect=T
model=setParam(model,'lb',{'R05807'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=4.99.1.3-RXN&redirect=T
model=setParam(model,'lb',{'R05808'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.1.1.151-RXN&redirect=T
model=setParam(model,'lb',{'R05809'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8761&redirect=T 
model=setParam(model,'lb',{'R05810'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8762&redirect=T
model=setParam(model,'lb',{'R05814'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8768&redirect=T
model=setParam(model,'lb',{'R05815'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8769&redirect=T
model=setParam(model,'lb',{'R05837'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13179
model=setParam(model,'lb',{'R05838'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=PDXJ-RXN&redirect=T
model=setParam(model,'lb',{'R05864'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12734
model=setParam(model,'lb',{'R05865'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3.7.1.9-RXN
model=setParam(model,'lb',{'R05884'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ISPH2-RXN&redirect=T
model=setParam(model,'lb',{'R06087'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3.2.1.48-RXN
model=setParam(model,'lb',{'R06088'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3.2.1.48-RXN 
model=setParam(model,'lb',{'R06129'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=TEICHOICSYN2-RXN
%R06171 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LTAA-RXN&redirect=T
                                    model=setParam(model,'lb',{'R06173'},[-1000]);model=setParam(model,'ub',{'R06173'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06173'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R06174'},[-1000]);model=setParam(model,'ub',{'R06174'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06174'))),:)=[3];  % Reaction directionality unclear!
%R06180 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TARTRATE-DEHYDROGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R06285'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.4.99.13&quickSearch=Quick+Search
model=setParam(model,'lb',{'R06366'},[0]);model=setParam(model,'ub',{'R06366'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14280&redirect=T
model=setParam(model,'lb',{'R06447'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8999&redirect=T
model=setParam(model,'lb',{'R06513'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DTDPGLUCDEHYDRAT-RXN&redirect=T
model=setParam(model,'lb',{'R06514'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DTDPDEHYDRHAMEPIM-RXN&redirect=T
model=setParam(model,'lb',{'R06529'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-6261&redirect=T
model=setParam(model,'lb',{'R06558'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14063&redirect=T
%R06590 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=2.2.1.1&quickSearch=Quick+Search
model=setParam(model,'lb',{'R06601'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3.5.2.17-RXN&redirect=T
model=setParam(model,'lb',{'R06633'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-15300&redirect=T 
model=setParam(model,'lb',{'R06786'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=MHPHYDROXY-RXN&redirect=T
model=setParam(model,'lb',{'R06787'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10040&redirect=T 
model=setParam(model,'lb',{'R06835'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9890&redirect=T
model=setParam(model,'lb',{'R06838'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9868&redirect=T 
model=setParam(model,'lb',{'R06895'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=HEMN-RXN&redirect=T
                                    model=setParam(model,'lb',{'R06897'},[-1000]);model=setParam(model,'ub',{'R06897'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06897'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R06909'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NAPHTHALENE-12-DIOXYGENASE-RXN
model=setParam(model,'lb',{'R06915'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.13.1&quickSearch=Quick+Search
                                    model=setParam(model,'lb',{'R06917'},[-1000]);model=setParam(model,'ub',{'R06917'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06917'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R06927'},[-1000]);model=setParam(model,'ub',{'R06927'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06927'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R06930'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.12.12&quickSearch=Quick+Search
model=setParam(model,'lb',{'R06936'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.13.1&quickSearch=Quick+Search
model=setParam(model,'lb',{'R06937'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.12.12&quickSearch=Quick+Search
model=setParam(model,'lb',{'R06939'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10078&redirect=T
%R06941 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-2044&redirect=T
%R06942 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-2425&redirect=T
model=setParam(model,'lb',{'R06983'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-2962
model=setParam(model,'lb',{'R07002'},[0]);model=setParam(model,'ub',{'R07002'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07003'},[0]);model=setParam(model,'ub',{'R07003'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07004'},[0]);model=setParam(model,'ub',{'R07004'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07023'},[0]);model=setParam(model,'ub',{'R07023'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07024'},[0]);model=setParam(model,'ub',{'R07024'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07025'},[0]);model=setParam(model,'ub',{'R07025'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07026'},[0]);model=setParam(model,'ub',{'R07026'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07034'},[0]);model=setParam(model,'ub',{'R07034'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=1.11.1.9&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07035'},[0]);model=setParam(model,'ub',{'R07035'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=1.11.1.9&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07069'},[0]);model=setParam(model,'ub',{'R07069'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07070'},[0]);model=setParam(model,'ub',{'R07070'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07083'},[0]);model=setParam(model,'ub',{'R07083'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07084'},[0]);model=setParam(model,'ub',{'R07084'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07091'},[0]);model=setParam(model,'ub',{'R07091'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07092'},[0]);model=setParam(model,'ub',{'R07092'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07093'},[0]);model=setParam(model,'ub',{'R07093'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07094'},[0]);model=setParam(model,'ub',{'R07094'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07100'},[0]);model=setParam(model,'ub',{'R07100'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
                                    model=setParam(model,'lb',{'R07105'},[-1000]);model=setParam(model,'ub',{'R07105'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07105'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R07113'},[-1000]);model=setParam(model,'ub',{'R07113'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07113'))),:)=[3];  % Reaction directionality unclear!  
                                    model=setParam(model,'lb',{'R07116'},[-1000]);model=setParam(model,'ub',{'R07116'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07116'))),:)=[3];  % Reaction directionality unclear!                                     
%R07136 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=R230-RXN&redirect=T
model=changeRxns(model,'R07168','NADH[c] + H+[c] + 5,10-Methylenetetrahydrofolate[c] <=> NAD+[c] + 5-Methyltetrahydrofolate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                           
model=setParam(model,'lb',{'R07168'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.5.1.20-RXN&redirect=T
model=setParam(model,'lb',{'R07210'},[0]);model=setParam(model,'ub',{'R07210'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-280&redirect=T
model=setParam(model,'lb',{'R07211'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN66-569&redirect=T
model=setParam(model,'lb',{'R07262'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9311&redirect=T
model=setParam(model,'lb',{'R07263'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=NAPHTHOATE-SYN-RXN&redirect=T
model=setParam(model,'lb',{'R07268'},[0]);model=setParam(model,'ub',{'R07268'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=BTUR2-RXN&redirect=T
model=setParam(model,'lb',{'R07270'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXNARA-8002&redirect=T
model=setParam(model,'lb',{'R07274'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2.5.1.65-RXN&redirect=T
model=setParam(model,'lb',{'R07280'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RIBOPHOSPHAT-RXN&redirect=T
model=setParam(model,'lb',{'R07281'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DIOHBUTANONEPSYN-RXN&redirect=T
model=setParam(model,'lb',{'R07302'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=6.3.1.10-RXN&redirect=T
%R07396 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=2.6.1.57&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07404'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-742&redirect=T
model=setParam(model,'lb',{'R07405'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-743&redirect=T
model=setParam(model,'lb',{'R07411'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=HEMEOSYN-RXN&redirect=T
                                    model=setParam(model,'lb',{'R07412'},[-1000]);model=setParam(model,'ub',{'R07412'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07412'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R07415'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-3921&redirect=T
                                    model=setParam(model,'lb',{'R07443'},[-1000]);model=setParam(model,'ub',{'R07443'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07443'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R07459'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9789&redirect=T
model=setParam(model,'lb',{'R07460'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-308&redirect=T
                                    model=setParam(model,'lb',{'R07463'},[-1000]);model=setParam(model,'ub',{'R07463'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07463'))),:)=[3];  % Reaction directionality unclear!
model=changeRxns(model,'R07605','2 NADPH[c] + 2 H+[c] + 7-Cyano-7-carbaguanine[c] <=> 2 NADP+[c] + 7-Aminomethyl-7-carbaguanine[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                               
model=setParam(model,'lb',{'R07605'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-4022&redirect=T
%R07618 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=1.8.1.4&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07644'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ENTMULTI-RXN&redirect=T
model=setParam(model,'lb',{'R07669'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.8.1.5&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07670'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.8.1.5&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07671'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-2461
model=setParam(model,'lb',{'R07704'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NAPHTHALENE-12-DIOXYGENASE-RXN
model=setParam(model,'lb',{'R07706'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-3661&redirect=T
model=setParam(model,'lb',{'R07709'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-10446
model=setParam(model,'lb',{'R07710'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-10446
model=setParam(model,'lb',{'R07762'},[0]);model=setParam(model,'ub',{'R07762'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.3.1.41&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07763'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.100&quickSearch=Quick+Search
                                    model=setParam(model,'lb',{'R07764'},[-1000]);model=setParam(model,'ub',{'R07764'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07764'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R07765'},[-1000]);model=setParam(model,'ub',{'R07765'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07765'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R07766'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-947&redirect=T
model=setParam(model,'lb',{'R07767'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-949&redirect=T
model=setParam(model,'lb',{'R07768'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.8.1.8&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07769'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.3.1.181&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07772'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8763&redirect=T
model=setParam(model,'lb',{'R07773'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8764&redirect=T
model=setParam(model,'lb',{'R07795'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13704
model=changeRxns(model,'R07891','CoA[c] + 3-Oxo-OPC8-CoA[c] <=> Acetyl-CoA[c] + OPC6-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                               
model=setParam(model,'lb',{'R07891'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.3.1.16&quickSearch=Quick+Search
model=changeRxns(model,'R07895','CoA[c] + 3-Oxo-OPC6-CoA[c] <=> Acetyl-CoA[c] + OPC4-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                               
model=setParam(model,'lb',{'R07895'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.3.1.16&quickSearch=Quick+Search
model=changeRxns(model,'R07899','CoA[c] + 3-Oxo-OPC4-CoA[c] <=> Acetyl-CoA[c] + (+)-7-Isojasmonic acid CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                               
model=setParam(model,'lb',{'R07899'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.3.1.16&quickSearch=Quick+Search
                                    model=setParam(model,'lb',{'R08014'},[-1000]);model=setParam(model,'ub',{'R08014'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08014'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R08017'},[-1000]);model=setParam(model,'ub',{'R08017'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08017'))),:)=[3];  % Reaction directionality unclear!                                     
                                    model=setParam(model,'lb',{'R08034'},[-1000]);model=setParam(model,'ub',{'R08034'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08034'))),:)=[3];  % Reaction directionality unclear!                                     
                                    model=setParam(model,'lb',{'R08042'},[-1000]);model=setParam(model,'ub',{'R08042'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08042'))),:)=[3];  % Reaction directionality unclear!                                     
model=setParam(model,'lb',{'R08056'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=4.2.1.40&quickSearch=Quick+Search
                                    model=setParam(model,'lb',{'R08087'},[-1000]);model=setParam(model,'ub',{'R08087'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08087'))),:)=[3];  % Reaction directionality unclear!  
model=setParam(model,'lb',{'R08089'},[0]);model=setParam(model,'ub',{'R08089'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11915&redirect=T 
model=setParam(model,'lb',{'R08090'},[0]);model=setParam(model,'ub',{'R08090'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=4.1.3.26-RXN&redirect=T
model=setParam(model,'lb',{'R08091'},[0]);model=setParam(model,'ub',{'R08091'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11917&redirect=T
model=setParam(model,'lb',{'R08093'},[0]);model=setParam(model,'ub',{'R08093'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11919&redirect=T
                                    model=setParam(model,'lb',{'R08094'},[-1000]);model=setParam(model,'ub',{'R08094'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08094'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R08095'},[0]);model=setParam(model,'ub',{'R08095'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11921&redirect=T
                                    model=setParam(model,'lb',{'R08096'},[-1000]);model=setParam(model,'ub',{'R08096'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08096'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R08111'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.25&quickSearch=Quick+Search
model=setParam(model,'lb',{'R08112'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.25&quickSearch=Quick+Search
model=setParam(model,'lb',{'R08113'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.25&quickSearch=Quick+Search
model=setParam(model,'lb',{'R08120'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.1.1.45&quickSearch=Quick+Search
model=setParam(model,'lb',{'R08121'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.1.1.45&quickSearch=Quick+Search
model=setParam(model,'lb',{'R08146'},[0]);model=setParam(model,'ub',{'R08146'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=1.2.1.3&quickSearch=Quick+Search
model=changeRxns(model,'R08210','2 H+[c] + 2 Reduced ferredoxin[c] + 1-Hydroxy-2-methyl-2-butenyl 4-diphosphate[c] <=> H2O[c] + 2 Oxidized ferredoxin[c] + Dimethylallyl diphosphate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                               
model=setParam(model,'lb',{'R08210'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.17.7.4&quickSearch=Quick+Search
model=setParam(model,'lb',{'R08218'},[0]);model=setParam(model,'ub',{'R08218'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-2161&redirect=T
model=setParam(model,'lb',{'R08222'},[-1000]);model=setParam(model,'ub',{'R08222'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THYM-PHOSPH-RXN                                         Irreversible in the original FDR.
%R08230 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=THYM-PHOSPH-RXN
model=setParam(model,'lb',{'R08240'},[-1000]);model=setParam(model,'ub',{'R08240'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=IMP-DEHYDROG-RXN                                        Irreversible in the original FDR.
model=setParam(model,'lb',{'R08244'},[0]);model=setParam(model,'ub',{'R08244'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GMP-SYN-GLUT-RXN
                                    model=setParam(model,'lb',{'R08280'},[-1000]);model=setParam(model,'ub',{'R08280'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08280'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R08281'},[-1000]);model=setParam(model,'ub',{'R08281'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08281'))),:)=[3];  % Reaction directionality unclear!    
                                    model=setParam(model,'lb',{'R08306'},[-1000]);model=setParam(model,'ub',{'R08306'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08306'))),:)=[3];  % Reaction directionality unclear!                                      
                                    model=setParam(model,'lb',{'R08310'},[-1000]);model=setParam(model,'ub',{'R08310'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08310'))),:)=[3];  % Reaction directionality unclear!                                      
model=setParam(model,'lb',{'R08359'},[0]);model=setParam(model,'ub',{'R08359'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-5217&redirect=T
%R08503 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14150&redirect=T
model=setParam(model,'lb',{'R08549'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=2OXOGLUTARATEDEH-RXN&redirect=T
                                    model=setParam(model,'lb',{'R08571'},[-1000]);model=setParam(model,'ub',{'R08571'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08571'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R08572'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GKI-RXN&redirect=T
model=setParam(model,'lb',{'R08618'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=R15-RXN&redirect=T
model=setParam(model,'lb',{'R08620'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-2203&redirect=T
model=setParam(model,'lb',{'R08624'},[0]);model=setParam(model,'ub',{'R08624'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXNQT-4164&redirect=T
model=setParam(model,'lb',{'R08628'},[0]);model=setParam(model,'ub',{'R08628'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXNQT-4167&redirect=T
model=setParam(model,'lb',{'R08634'},[0]);model=setParam(model,'ub',{'R08634'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXNQT-4170&redirect=T
%R08639 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHOSPHOGLUCMUT-RXN
model=setParam(model,'lb',{'R08641'},[0]);model=setParam(model,'ub',{'R08641'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXNQT-4173&redirect=T
model=setParam(model,'lb',{'R08645'},[0]);model=setParam(model,'ub',{'R08645'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXNQT-4176&redirect=T
model=setParam(model,'lb',{'R08648'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=ACETOOHBUTSYN-RXN&redirect=T
model=setParam(model,'lb',{'R08689'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-882&redirect=T
model=setParam(model,'lb',{'R08702'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.5.3.-&quickSearch=Quick+Search
%R08714 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8&redirect=T
%R08734 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9797&redirect=T
%R08739 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9845&redirect=T
model=setParam(model,'lb',{'R08836'},[0]);model=setParam(model,'ub',{'R08836'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8460&redirect=T
model=setParam(model,'lb',{'R08871'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10837&redirect=T
model=setParam(model,'lb',{'R08938'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9382&redirect=T
model=setParam(model,'lb',{'R08968'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18662
%R09097 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12736&redirect=T
%R09099 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLYOHMETRANS-RXN
model=setParam(model,'lb',{'R09125'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1.13.11.9-RXN&redirect=T
model=setParam(model,'lb',{'R09136'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.1.1.45&quickSearch=Quick+Search
model=setParam(model,'lb',{'R09159'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NAPHTHALENE-12-DIOXYGENASE-RXN
model=setParam(model,'lb',{'R09220'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10107&redirect=T
model=setParam(model,'lb',{'R09222'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.1.1.45&quickSearch=Quick+Search
model=setParam(model,'lb',{'R09233'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8820&redirect=T
model=setParam(model,'lb',{'R09248'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8992&redirect=T
model=setParam(model,'lb',{'R09272'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10889&redirect=T
model=setParam(model,'lb',{'R09365'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12730&redirect=T
                                    model=setParam(model,'lb',{'R09372'},[-1000]);model=setParam(model,'ub',{'R09372'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09372'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R09381'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=1-ACYLGLYCEROL-3-P-ACYLTRANSFER-RXN&redirect=T
model=setParam(model,'lb',{'R09382'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=TRNA-CYTIDYLYLTRANSFERASE-RXN&redirect=T
                                    model=setParam(model,'lb',{'R09383'},[-1000]);model=setParam(model,'ub',{'R09383'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09383'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R09384'},[-1000]);model=setParam(model,'ub',{'R09384'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09384'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R09386'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.7.7.72&quickSearch=Quick+Search
model=setParam(model,'lb',{'R09394'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8340&redirect=T
                                    model=setParam(model,'lb',{'R09395'},[-1000]);model=setParam(model,'ub',{'R09395'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09395'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R09409'},[0]);model=setParam(model,'ub',{'R09409'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.18&quickSearch=Quick+Search
model=setParam(model,'lb',{'R09493'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11193
model=setParam(model,'lb',{'R09497'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.7.5.1&quickSearch=Quick+Search 
model=setParam(model,'lb',{'R09513'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9770&redirect=T
model=setParam(model,'lb',{'R09518'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12219
model=setParam(model,'lb',{'R09543'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11475&redirect=T
model=setParam(model,'lb',{'R09554'},[0]);model=setParam(model,'ub',{'R09554'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11053&redirect=T
model=setParam(model,'lb',{'R09555'},[0]);model=setParam(model,'ub',{'R09555'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-2401&redirect=T
model=setParam(model,'lb',{'R09556'},[0]);model=setParam(model,'ub',{'R09556'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11054&redirect=T
model=setParam(model,'lb',{'R09565'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GALLATE-DIOXYGENASE-RXN
model=setParam(model,'lb',{'R09597'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-1961&redirect=T
model=setParam(model,'lb',{'R09726'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8344
model=setParam(model,'lb',{'R09735'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8348&redirect=T
model=setParam(model,'lb',{'R09736'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11046&redirect=T
model=setParam(model,'lb',{'R09763'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.4.99.14&quickSearch=Quick+Search
model=setParam(model,'lb',{'R09764'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.4.99.14&quickSearch=Quick+Search
model=setParam(model,'lb',{'R09766'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.4.99.15&quickSearch=Quick+Search
model=setParam(model,'lb',{'R09768'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=5.3.1.28&quickSearch=Quick+Search
model=setParam(model,'lb',{'R09769'},[0]);model=setParam(model,'ub',{'R09769'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=5.3.1.28&quickSearch=Quick+Search
model=setParam(model,'lb',{'R09771'},[0]);model=setParam(model,'ub',{'R09771'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11365&redirect=T
%R09837 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-6510
model=setParam(model,'lb',{'R09838'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-2042
model=setParam(model,'lb',{'R09839'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-6512
                                    model=setParam(model,'lb',{'R09840'},[-1000]);model=setParam(model,'ub',{'R09840'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09840'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R09845'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11698
model=setParam(model,'lb',{'R09978'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12093
model=setParam(model,'lb',{'R10035'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.2.1.21&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10039'},[0]);model=setParam(model,'ub',{'R10039'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=3.2.1.21&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10040'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.2.1.21&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10042'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17503
model=setParam(model,'lb',{'R10043'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17503
model=setParam(model,'lb',{'R10052'},[0]);model=setParam(model,'ub',{'R10052'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.85&quickSearch=Quick+Search
%R10074 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=GLUTARYL-COA-DEHYDROG-RXN
%R10092 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-5224
model=setParam(model,'lb',{'R10115'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11474
model=setParam(model,'lb',{'R10116'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.100&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10117'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=4.2.1.59&quickSearch=Quick+Search 
model=setParam(model,'lb',{'R10118'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.10&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10119'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.3.1.41&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10120'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.1.1.100&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10121'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=4.2.1.59&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10122'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.3.1.10&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10124'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11484
                                    model=setParam(model,'lb',{'R10125'},[-1000]);model=setParam(model,'ub',{'R10125'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10125'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R10126'},[-1000]);model=setParam(model,'ub',{'R10126'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10126'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R10147'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DIHYDRODIPICSYN-RXN&redirect=T
                                    model=setParam(model,'lb',{'R10151'},[-1000]);model=setParam(model,'ub',{'R10151'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10151'))),:)=[3];  % Reaction directionality unclear!
model=changeRxns(model,'R10159','NADPH[c] + 2 Oxidized adrenal ferredoxin[c] <=> NADP+[c] + H+[c] + 2 Reduced adrenal ferredoxin[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R10159'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-13685&redirect=T
model=setParam(model,'lb',{'R10170'},[-1000]);model=setParam(model,'ub',{'R10170'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13163                                           Irreversible in the original FDR.
model=setParam(model,'lb',{'R10177'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-13323&redirect=T
model=setParam(model,'lb',{'R10206'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18086
model=setParam(model,'lb',{'R10209'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-1321&redirect=T
                                    model=setParam(model,'lb',{'R10210'},[-1000]);model=setParam(model,'ub',{'R10210'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10210'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R10213'},[-1000]);model=setParam(model,'ub',{'R10213'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10213'))),:)=[3];  % Reaction directionality unclear!
%R10223 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-13997&redirect=T
model=setParam(model,'lb',{'R10247'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=THIAZOLSYN2-RXN&redirect=T
%R10343 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8807&redirect=T
%R10463 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14569&redirect=T
%R10619 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=5.1.3.3&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10645'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.8.4.3&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10646'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.8.4.3&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10647'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.8.4.3&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10648'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14570
model=setParam(model,'lb',{'R10652'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-6366
model=setParam(model,'lb',{'R10657'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14706
model=setParam(model,'lb',{'R10658'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14705
model=setParam(model,'lb',{'R10705'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17572
model=setParam(model,'lb',{'R10707'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=2.3.1.180-RXN
model=setParam(model,'lb',{'R10712'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.3&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10757'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.5.1.74&quickSearch=Quick+Search
%R10764 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-5298
model=setParam(model,'lb',{'R10806'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14992
model=setParam(model,'lb',{'R10816'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.6.1.-&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10820'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12571
model=setParam(model,'lb',{'R10831'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.2.1.52&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10841'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=URONATE-DEHYDROGENASE-RXN
%R10845 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-14469
model=changeRxns(model,'R10859','Reduced flavodoxin[c] + 2-C-Methyl-D-erythritol 2,4-cyclodiphosphate[c] <=> H2O[c] + Oxidized flavodoxin[c] + 1-Hydroxy-2-methyl-2-butenyl 4-diphosphate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R10859'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15878
model=setParam(model,'lb',{'R10907'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16652
model=setParam(model,'lb',{'R10936'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12573
                                    model=setParam(model,'lb',{'R10948'},[-1000]);model=setParam(model,'ub',{'R10948'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10948'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R10949'},[-1000]);model=setParam(model,'ub',{'R10949'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10949'))),:)=[3];  % Reaction directionality unclear!
%R10991 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=2.6.1.42&quickSearch=Quick+Search
model=setParam(model,'lb',{'R10993'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-20400
model=setParam(model,'lb',{'R10994'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-20401
model=setParam(model,'lb',{'R11018'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-7563
model=setParam(model,'lb',{'R11024'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18660
model=setParam(model,'lb',{'R11025'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18661
model=setParam(model,'lb',{'R11037'},[-1000]);model=setParam(model,'ub',{'R11037'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16480                                   Irreversible in the original FDR.
                                    model=setParam(model,'lb',{'R11062'},[-1000]);model=setParam(model,'ub',{'R11062'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R11062'))),:)=[3];  % Reaction directionality unclear!
%R11073 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-10856&redirect=T
model=setParam(model,'lb',{'R11130'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PROPCOASYN-RXN
model=setParam(model,'lb',{'R11168'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17357
model=setParam(model,'lb',{'R11173'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=3.1.3.73&quickSearch=Quick+Search
model=setParam(model,'lb',{'R11174'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.7.8.26&quickSearch=Quick+Search
model=setParam(model,'lb',{'R11188'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16009
model=setParam(model,'lb',{'R11225'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16937
%R11319 is reversible                                                                       %https://biocyc.org/META/substring-search?type=NIL&object=2.7.4.3&quickSearch=Quick+Search
model=changeRxns(model,'R11329','Coproporphyrin III[c] + Fe2+[c] <=> 2 H+[c] + Fe-coproporphyrin III[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R11329'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17518
model=setParam(model,'lb',{'R11372'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17809
model=setParam(model,'lb',{'R11396'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16748
model=setParam(model,'lb',{'R11443'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=LEUCYLTRANSFERASE-RXN
model=setParam(model,'lb',{'R11444'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17846
model=setParam(model,'lb',{'R11528'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.8.1.7&quickSearch=Quick+Search
%R11529 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12587&redirect=T
model=setParam(model,'lb',{'R11546'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17813
model=setParam(model,'lb',{'R11547'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17898
model=setParam(model,'lb',{'R11548'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17899
model=setParam(model,'lb',{'R11555'},[0]);model=setParam(model,'ub',{'R11555'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16283
model=setParam(model,'lb',{'R11556'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16282
model=setParam(model,'lb',{'R11557'},[0]);model=setParam(model,'ub',{'R11557'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16284
model=setParam(model,'lb',{'R11581'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-262
model=setParam(model,'lb',{'R11582'},[0]);model=setParam(model,'ub',{'R11582'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-6254
model=setParam(model,'lb',{'R11633'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-746
model=setParam(model,'lb',{'R11634'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-745
model=setParam(model,'lb',{'R11635'},[0]);model=setParam(model,'ub',{'R11635'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-724
model=setParam(model,'lb',{'R11636'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-723
                                    model=setParam(model,'lb',{'R11671'},[-1000]);model=setParam(model,'ub',{'R11671'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R11671'))),:)=[3];  % Reaction directionality unclear!
model=changeRxns(model,'R11765','NADPH[c] + H+[c] + 7,8-Dihydrobiopterin[c] <=> NADP+[c] + Tetrahydrobiopterin[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R11765'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-19329
model=setParam(model,'lb',{'R11768'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8775
model=setParam(model,'lb',{'R11785'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18659
model=setParam(model,'lb',{'R11786'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15952
model=setParam(model,'lb',{'R11861'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-19024
model=setParam(model,'lb',{'R11894'},[0]);model=setParam(model,'ub',{'R11894'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=UDPKIN-RXN
model=setParam(model,'lb',{'R11895'},[0]);model=setParam(model,'ub',{'R11895'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DUDPKIN-RXN
model=setParam(model,'lb',{'R11896'},[0]);model=setParam(model,'ub',{'R11896'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=DUTP-PYROP-RXN
model=setParam(model,'lb',{'R11906'},[0]);model=setParam(model,'ub',{'R11906'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8667
model=setParam(model,'lb',{'R12024'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18696
model=setParam(model,'lb',{'R12049'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15021
model=setParam(model,'lb',{'R12096'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11952
model=setParam(model,'lb',{'R12097'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18797
model=setParam(model,'lb',{'R12173'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18797
model=setParam(model,'lb',{'R12188'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-6763
                                    model=setParam(model,'lb',{'R12193'},[-1000]);model=setParam(model,'ub',{'R12193'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R12193'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R12202'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.7.8.43&quickSearch=Quick+Search
%R12240 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9590
model=setParam(model,'lb',{'R12241'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9591
model=setParam(model,'lb',{'R12423'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.8.1.8&quickSearch=Quick+Search
model=setParam(model,'lb',{'R12424'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=2.8.1.8&quickSearch=Quick+Search
model=setParam(model,'lb',{'R12427'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-947
model=setParam(model,'lb',{'R12428'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13037
model=changeRxns(model,'R00647','L-xylo-Hexulonolactone[c] <=> Ascorbate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R00647'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8784&redirect=T
model=changeRxns(model,'R01652','(2S)-2-Isopropyl-3-oxosuccinate[c] <=> CO2[c] + 4-Methyl-2-oxopentanoate[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R01652'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-7800&redirect=T
%R02317 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8173&redirect=T
model=setParam(model,'lb',{'R02962'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8483&redirect=T
                                    model=setParam(model,'lb',{'R03166'},[-1000]);model=setParam(model,'ub',{'R03166'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R03166'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R03186'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-3523&redirect=T
model=changeRxns(model,'R03672','Dopaquinone[c] + 2-Carboxy-2,3-dihydro-5,6-dihydroxyindole[c] <=> 3,4-Dihydroxy-L-phenylalanine[c] + L-Dopachrome[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R03672'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11369&redirect=T
model=setParam(model,'lb',{'R03674'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-11403&redirect=T
model=changeRxns(model,'R03758','L-2-Amino-3-oxobutanoic acid[c] <=> CO2[c] + Aminoacetone[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R03758'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=THREOSPON-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04500'},[-1000]);model=setParam(model,'ub',{'R04500'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04500'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R04578'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=C04718-C060066-AUTOTRANSFORMATION-RXN&redirect=T
                                    model=setParam(model,'lb',{'R04639'},[-1000]);model=setParam(model,'ub',{'R04639'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04639'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R04701'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8181&redirect=T
model=setParam(model,'lb',{'R04769'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8144&redirect=T
                                    model=setParam(model,'lb',{'R05041'},[-1000]);model=setParam(model,'ub',{'R05041'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05041'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R05048'},[-1000]);model=setParam(model,'ub',{'R05048'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05048'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R05098'},[-1000]);model=setParam(model,'ub',{'R05098'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05098'))),:)=[3];  % Reaction directionality unclear!  
                                    model=setParam(model,'lb',{'R05455'},[-1000]);model=setParam(model,'ub',{'R05455'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05455'))),:)=[3];  % Reaction directionality unclear!                                    
model=changeRxns(model,'R05483','1,3,4,6-Tetrachloro-1,4-cyclohexadiene[c] <=> Hydrochloric acid[c] + 1,2,4-Trichlorobenzene[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R05483'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LINABSPON-RXN&redirect=T
model=setParam(model,'lb',{'R05484'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=LINBBSPON-RXN&redirect=T
                                    model=setParam(model,'lb',{'R05824'},[-1000]);model=setParam(model,'ub',{'R05824'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05824'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R05826'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8225&redirect=T
                                    model=setParam(model,'lb',{'R05844'},[-1000]);model=setParam(model,'ub',{'R05844'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05844'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R05845'},[-1000]);model=setParam(model,'ub',{'R05845'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05845'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R05846'},[-1000]);model=setParam(model,'ub',{'R05846'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05846'))),:)=[3];  % Reaction directionality unclear!                                     
                                    model=setParam(model,'lb',{'R05887'},[-1000]);model=setParam(model,'ub',{'R05887'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05887'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R05889'},[-1000]);model=setParam(model,'ub',{'R05889'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05889'))),:)=[3];  % Reaction directionality unclear!   
                                    model=setParam(model,'lb',{'R06063'},[-1000]);model=setParam(model,'ub',{'R06063'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06063'))),:)=[3];  % Reaction directionality unclear!                                     
                                    model=setParam(model,'lb',{'R06064'},[-1000]);model=setParam(model,'ub',{'R06064'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06064'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R06194'},[-1000]);model=setParam(model,'ub',{'R06194'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06194'))),:)=[3];  % Reaction directionality unclear!                                     
                                    model=setParam(model,'lb',{'R06650'},[-1000]);model=setParam(model,'ub',{'R06650'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06650'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R06830'},[-1000]);model=setParam(model,'ub',{'R06830'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06830'))),:)=[3];  % Reaction directionality unclear!                                     
model=setParam(model,'lb',{'R06982'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-2961&redirect=T
                                    model=setParam(model,'lb',{'R07009'},[-1000]);model=setParam(model,'ub',{'R07009'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07009'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R07010'},[-1000]);model=setParam(model,'ub',{'R07010'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07010'))),:)=[3];  % Reaction directionality unclear!                                    
                                    model=setParam(model,'lb',{'R07011'},[-1000]);model=setParam(model,'ub',{'R07011'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07011'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R07012'},[-1000]);model=setParam(model,'ub',{'R07012'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07012'))),:)=[3];  % Reaction directionality unclear!
model=setParam(model,'lb',{'R07198'},[0]);                                                  %https://biocyc.org/META/substring-search?type=NIL&object=1.14.14.162&quickSearch=Quick+Search
model=setParam(model,'lb',{'R07316'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN0-5222&redirect=T
model=setParam(model,'lb',{'R07406'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13179
%R07408 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-6423&redirect=T
model=setParam(model,'lb',{'R07679'},[0]);model=setParam(model,'ub',{'R07679'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12876&redirect=T
                                    model=setParam(model,'lb',{'R07947'},[-1000]);model=setParam(model,'ub',{'R07947'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07947'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R07948'},[-1000]);model=setParam(model,'ub',{'R07948'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R07948'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R08020'},[-1000]);model=setParam(model,'ub',{'R08020'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08020'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R08279'},[-1000]);model=setParam(model,'ub',{'R08279'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08279'))),:)=[3];  % Reaction directionality unclear!                                    
                                    model=setParam(model,'lb',{'R08289'},[-1000]);model=setParam(model,'ub',{'R08289'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08289'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R08358'},[-1000]);model=setParam(model,'ub',{'R08358'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08358'))),:)=[3];  % Reaction directionality unclear!   
                                    model=setParam(model,'lb',{'R08361'},[-1000]);model=setParam(model,'ub',{'R08361'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08361'))),:)=[3];  % Reaction directionality unclear!                                    
                                    model=setParam(model,'lb',{'R08372'},[-1000]);model=setParam(model,'ub',{'R08372'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08372'))),:)=[3];  % Reaction directionality unclear!
                                    model=setParam(model,'lb',{'R08637'},[-1000]);model=setParam(model,'ub',{'R08637'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08637'))),:)=[3];  % Reaction directionality unclear!                                    
model=setParam(model,'lb',{'R08698'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8899&redirect=T
model=setParam(model,'lb',{'R08818'},[0]);model=setParam(model,'ub',{'R08818'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8485&redirect=T
model=setParam(model,'lb',{'R08820'},[0]);model=setParam(model,'ub',{'R08820'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8482&redirect=T
model=setParam(model,'lb',{'R08821'},[0]);model=setParam(model,'ub',{'R08821'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8478&redirect=T
model=setParam(model,'lb',{'R08825'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8487&redirect=T
model=setParam(model,'lb',{'R08827'},[0]);model=setParam(model,'ub',{'R08827'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8491&redirect=T
model=setParam(model,'lb',{'R08829'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8490&redirect=T
model=setParam(model,'lb',{'R08831'},[0]);model=setParam(model,'ub',{'R08831'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8567&redirect=T
model=setParam(model,'lb',{'R08833'},[0]);model=setParam(model,'ub',{'R08833'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8624&redirect=T
model=setParam(model,'lb',{'R08834'},[0]);model=setParam(model,'ub',{'R08834'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8625&redirect=T
model=setParam(model,'lb',{'R08835'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8486&redirect=T
model=setParam(model,'lb',{'R08837'},[0]);model=setParam(model,'ub',{'R08837'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-8461&redirect=T
model=setParam(model,'lb',{'R08976'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-13939&redirect=T
                                    model=setParam(model,'lb',{'R09142'},[-1000]);model=setParam(model,'ub',{'R09142'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09142'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R09143'},[-1000]);model=setParam(model,'ub',{'R09143'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09143'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R09144'},[-1000]);model=setParam(model,'ub',{'R09144'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09144'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R09145'},[-1000]);model=setParam(model,'ub',{'R09145'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09145'))),:)=[3];  % Reaction directionality unclear! 
%R09149 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=MOXXANAU-RXN
                                    model=setParam(model,'lb',{'R09274'},[-1000]);model=setParam(model,'ub',{'R09274'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09274'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R09367'},[0]);model=setParam(model,'ub',{'R09367'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12864
                                    model=setParam(model,'lb',{'R09368'},[-1000]);model=setParam(model,'ub',{'R09368'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09368'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R09391'},[0]);model=setParam(model,'ub',{'R09391'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-9418&redirect=T
model=setParam(model,'lb',{'R09392'},[0]);model=setParam(model,'ub',{'R09392'},[1000]);     %https://biocyc.org/META/substring-search?type=NIL&object=1.14.13.105&quickSearch=Quick+Search
                                    model=setParam(model,'lb',{'R09431'},[-1000]);model=setParam(model,'ub',{'R09431'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09431'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R09432'},[-1000]);model=setParam(model,'ub',{'R09432'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09432'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R09433'},[-1000]);model=setParam(model,'ub',{'R09433'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09433'))),:)=[3];  % Reaction directionality unclear! 
                                    model=setParam(model,'lb',{'R09434'},[-1000]);model=setParam(model,'ub',{'R09434'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R09434'))),:)=[3];  % Reaction directionality unclear!                                     
model=setParam(model,'lb',{'R09885'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12715
model=setParam(model,'lb',{'R09981'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN0-6450
model=setParam(model,'lb',{'R10013'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13085
                                    model=setParam(model,'lb',{'R10018'},[-1000]);model=setParam(model,'ub',{'R10018'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10018'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R10022'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13044
model=setParam(model,'lb',{'R10030'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9640
model=setParam(model,'lb',{'R10034'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-5066
model=setParam(model,'lb',{'R10041'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12065
model=setParam(model,'lb',{'R10078'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12799
model=setParam(model,'lb',{'R10082'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-12536&redirect=T
model=setParam(model,'lb',{'R10097'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16227
model=setParam(model,'lb',{'R10245'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=1.4.3.19-RXN
model=setParam(model,'lb',{'R10302'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-6581
model=setParam(model,'lb',{'R10792'},[0]);model=setParam(model,'ub',{'R10792'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15673
model=setParam(model,'lb',{'R10854'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15852
model=setParam(model,'lb',{'R10856'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15853
model=setParam(model,'lb',{'R10914'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16319
model=setParam(model,'lb',{'R11017'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15850
model=setParam(model,'lb',{'R11034'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12943
model=setParam(model,'lb',{'R11049'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16527
model=setParam(model,'lb',{'R11098'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15123
model=setParam(model,'lb',{'R11099'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15127
model=setParam(model,'lb',{'R11100'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15124
model=setParam(model,'lb',{'R11101'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15121
model=changeRxns(model,'R11148','3,4-Dihydroxy-2-methy-4-[(2E,6E)-farnesyl]-3,4-dihydroquinoline 1-oxide[c] <=> H2O[c] + Aurachin B[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
model=setParam(model,'lb',{'R11148'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17336
model=setParam(model,'lb',{'R11176'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-16527
model=setParam(model,'lb',{'R11228'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17253
                                    model=setParam(model,'lb',{'R11523'},[-1000]);model=setParam(model,'ub',{'R11523'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R11523'))),:)=[3];  % Reaction directionality unclear! 
model=setParam(model,'lb',{'R11550'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12350
model=setParam(model,'lb',{'R11570'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17689
%R11571 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17692
model=setParam(model,'lb',{'R11572'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17695
model=setParam(model,'lb',{'R11573'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17696
%R11574 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12020
model=setParam(model,'lb',{'R11575'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17687
model=setParam(model,'lb',{'R11576'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17691
model=setParam(model,'lb',{'R11577'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17694
model=setParam(model,'lb',{'R11691'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13809
model=setParam(model,'lb',{'R11695'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13811
model=setParam(model,'lb',{'R11709'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-17455
model=setParam(model,'lb',{'R11736'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18158
model=setParam(model,'lb',{'R11739'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-11420
model=setParam(model,'lb',{'R11865'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-18089
model=setParam(model,'lb',{'R11947'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13161
model=setParam(model,'lb',{'R12091'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13477
model=setParam(model,'lb',{'R12094'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-13469
model=setParam(model,'lb',{'R12132'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-8971
                                    model=setParam(model,'lb',{'R12140'},[-1000]);model=setParam(model,'ub',{'R12140'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R12140'))),:)=[3];  % Reaction directionality unclear! 
%R12185 is reversible                                                                       %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-19659
model=setParam(model,'lb',{'R12268'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-20054
model=setParam(model,'lb',{'R12274'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-20050
model=setParam(model,'lb',{'R12310'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-12965
model=setParam(model,'lb',{'R12335'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-20066
model=setParam(model,'lb',{'R12421'},[0]);                                                  %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15750


%% Manual curation of reaction directionalities pertaining to reactions added in conjunction with metabolites in category 2

model=setParam(model,'lb',{'R00321'},[0]);model=setParam(model,'ub',{'R00321'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-14728&redirect=T
model=setParam(model,'lb',{'R02526'},[0]);model=setParam(model,'ub',{'R02526'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=L-ARABINONOLACTONASE-RXN&redirect=T
                                    model=setParam(model,'lb',{'R01757'},[-1000]);model=setParam(model,'ub',{'R01757'},[1000]);                                                                      % Reaction directionality unclear! 
model=setParam(model,'lb',{'R01176'},[0]);model=setParam(model,'ub',{'R01176'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-19193
model=setParam(model,'lb',{'R03544'},[-1000]);model=setParam(model,'ub',{'R03544'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-161&redirect=T
model=setParam(model,'lb',{'R01961'},[0]);model=setParam(model,'ub',{'R01961'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=GLUCOSAMINE-KINASE-RXN&redirect=T
model=setParam(model,'lb',{'R01628'},[0]);model=setParam(model,'ub',{'R01628'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3-HYDROXYBENZOATE-4-MONOOXYGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R02589'},[0]);model=setParam(model,'ub',{'R02589'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=3-HYDROXYBENZOATE-6-MONOOXYGENASE-RXN&redirect=T
model=setParam(model,'lb',{'R01100'},[0]);model=setParam(model,'ub',{'R01100'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=BETAGALACTOSID-RXN
                                    model=setParam(model,'lb',{'R01898'},[-1000]);model=setParam(model,'ub',{'R01898'},[1000]);                                                                      % Reaction directionality unclear! 
model=setParam(model,'lb',{'R00877'},[-1000]);model=setParam(model,'ub',{'R00877'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=MANNOSE-ISOMERASE-RXN&redirect=T
model=setParam(model,'lb',{'R02377'},[-1000]);model=setParam(model,'ub',{'R02377'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=13-PROPANEDIOL-DEHYDROGENASE-RXN
model=setParam(model,'lb',{'R00875'},[-1000]);model=setParam(model,'ub',{'R00875'},[1000]); %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=RXN-7644&redirect=T
model=setParam(model,'lb',{'R00010'},[0]);model=setParam(model,'ub',{'R00010'},[1000]);     %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=TREHALA-RXN

%% Manual curation of reaction directionalities pertaining to reactions labelled as 'incomplete' or reactions w/ 'undefined stoichiometries'

                                                                            % NB! The following reactions appears in the FDR upon generating a model with the command
                                                                            % "model=getKEGGModelForOrganism('hpse');" which includes reactions with undefined stoichiometries 
                                                                            % etc. As such, they are kept here in case users prefer generating an FDR using these settings. 
                                                                            % If this is not the case, the following constraints are to be kept silenced.

                                    %model=setParam(model,'lb',{'R08763'},[-1000]);model=setParam(model,'ub',{'R08763'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08763'))),:)=[3];  % Reaction directionality unclear!
%model=setParam(model,'lb',{'R00381'},[0]);model=setParam(model,'ub',{'R00381'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00381'))),:)=[4];    %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DNA-LIGASE-ATP-RXN&redirect=T
%model=setParam(model,'lb',{'R00382'},[0]);model=setParam(model,'ub',{'R00382'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00382'))),:)=[4];    %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=DNA-LIGASE-NAD%2b-RXN&redirect=T
%model=setParam(model,'lb',{'R02887'},[0]);model=setParam(model,'ub',{'R02887'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R02887'))),:)=[4];    %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=3.2.1.91-RXN
%model=setParam(model,'lb',{'R04241'},[0]);model=setParam(model,'ub',{'R04241'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04241'))),:)=[4];    %https://biocyc.org/META/NEW-IMAGE?type=NIL&object=FOLYLPOLYGLUTAMATESYNTH-RXN&redirect=T
%model=setParam(model,'lb',{'R04254'},[0]);model=setParam(model,'ub',{'R04254'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04254'))),:)=[4];    %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN1-42
                                    %model=setParam(model,'lb',{'R05196'},[-1000]);model=setParam(model,'ub',{'R05196'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R05196'))),:)=[3];  % Reaction directionality unclear! 
%model=setParam(model,'lb',{'R10152'},[0]);model=setParam(model,'ub',{'R10152'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10152'))),:)=[4];    %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=R17-RXN
%model=setParam(model,'lb',{'R10822'},[0]);model=setParam(model,'ub',{'R10822'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10822'))),:)=[4];    %https://biocyc.org/META/substring-search?type=NIL&object=6.5.1.7&quickSearch=Quick+Search
%model=setParam(model,'lb',{'R10823'},[0]);model=setParam(model,'ub',{'R10823'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R10823'))),:)=[4];    %https://biocyc.org/META/substring-search?type=NIL&object=6.5.1.7&quickSearch=Quick+Search
%model=setParam(model,'lb',{'R00590'},[0]);model=setParam(model,'ub',{'R00590'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00590'))),:)=[4];    %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-15125
%model=setParam(model,'lb',{'R00698'},[0]);model=setParam(model,'ub',{'R00698'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R00698'))),:)=[4];    %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=PHENYLALANINE-2-MONOOXYGENASE-RXN
%model=changeRxns(model,'R04546','CoA[c] + 3alpha,7alpha-Dihydroxy-5beta-cholestanoyl-CoA[c] <=> Propanoyl-CoA[c] + Chenodeoxycholoyl-CoA[c]',3);  % Reverses the order in which reactants and products appear in the reaction equation.                                                                                                                                                                                                                   
%model=setParam(model,'lb',{'R04546'},[0]);model=setParam(model,'lb',{'R04546'},[-1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R04546'))),:)=[4];   %https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=RXN-9850
%model=setParam(model,'lb',{'R06411'},[0]);model=setParam(model,'ub',{'R06411'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06411'))),:)=[4];    %https://biocyc.org/META/substring-search?type=NIL&object=4.2.1.17&quickSearch=Quick+Search
%model=setParam(model,'lb',{'R06412'},[0]);model=setParam(model,'ub',{'R06412'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R06412'))),:)=[4];    %https://biocyc.org/META/substring-search?type=NIL&object=4.2.1.17&quickSearch=Quick+Search
%model=setParam(model,'lb',{'R08700'},[0]);model=setParam(model,'ub',{'R08700'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08700'))),:)=[4];    %https://biocyc.org/META/substring-search?type=NIL&object=2.8.1.13&quickSearch=Quick+Search
                                    %model=setParam(model,'lb',{'R08701'},[-1000]);model=setParam(model,'ub',{'R08701'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R08701'))),:)=[3];  % Reaction directionality unclear!
                                    %model=setParam(model,'lb',{'R11905'},[-1000]);model=setParam(model,'ub',{'R11905'},[1000]);model.rxnConfidenceScores((find(ismember(model.rxns,'R11905'))),:)=[3];  % Reaction directionality unclear!                                     

                                                                            %% Artificial biomass reaction
                                   
                                                                            % The artificial biomass reaction adapted for the purpose of simulating production of 
                                                                            % biomass – i.e., growth – was adapted from the GEM RehMBEL1391_sbml_L3V1 on Cupriavidus 
                                                                            % necator H16 (formerly Ralstonia eutropha H16) originally reported by Park et al. (2011) 
                                                                            % and available in an updated version from the GitHub-repository of GitHub-user m-jahn 
                                                                            % (https://github.com/m-jahn/genome-scale-models). The adapted biomass reaction was designed 
                                                                            % as a lumping together of the following biomolecular pools: lipopolisaccharide, RNA, 
                                                                            % carbohydrate, phospholipid, peptidoglycan, protein, cofactors and vitamines and DNA. 
                                                                            % Accordingly, artificial biosynthesis reactions for each of these pools had to be adapted 
                                                                            % as well. All of these artificial reactions were assigned confidence score of 0 to reflect 
                                                                            % the fact that no computational or experimental efforts were done to check their eligibility 
                                                                            % in the specific case of H. pseudoflava. Metabolite IDs pertaining to biomass precursors 
                                                                            % were manually harmonized with KEGG-based terminology.
                                                                            

rxns.rxns={'lipopolisaccharide'};                                                       % Creates this new reaction ID.
rxns.equations={'0.28 CDP-ethanolamine + 0.14 Di[3-deoxy-D-manno-octulosonyl]-lipid IV(A) + 0.42 CMP-3-deoxy-D-manno-octulosonate + 0.42 ADP-L-glycero-D-manno-heptose + 0.28 UDP-glucose => lipopolisaccharide + 0.28 UDP + 0.42 CMP + 0.28 CDP + 0.42 ADP'};
model=addRxns(model,rxns,2,'c',true);                                                   % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'lipopolisaccharide'))),:)=[0];     % Assigns a confidence score of 0 as the reaction is artificial.


rxns.rxns={'RNA'};                                                                      % Creates this new reaction ID.
rxns.equations={'0.75 GTP + 0.747 UTP + 0.998 CTP + 0.631 ATP => RNA + 1.25 Orthophosphate + 1.25 ADP'};
model=addRxns(model,rxns,2,'c',true);                                                   % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'RNA'))),:)=[0];                    % Assigns a confidence score of 0 as the reaction is artificial.


rxns.rxns={'carbohydrate'};                                                             % Creates this new reaction ID.
rxns.equations={'0.984 UDP-N-acetyl-D-galactosamine + 3.937 UDP-N-acetyl-D-glucosamine => carbohydrate + 4.921 UDP'};
model=addRxns(model,rxns,2,'c',true);                                                   % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'carbohydrate'))),:)=[0];           % Assigns a confidence score of 0 as the reaction is artificial.


rxns.rxns={'phospholipid'};                                                             % Creates this new reaction ID.
rxns.equations={'0.927 Phosphatidylethanolamine + 0.283 Phosphatidylglycerol + 0.093 Cardiolipin => phospholipid'};
model=addRxns(model,rxns,2,'c',true);                                                   % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'phospholipid'))),:)=[0];           % Assigns a confidence score of 0 as the reaction is artificial.


rxns.rxns={'peptidoglycan'};                                                           % Creates this new reaction ID.
rxns.equations={'D-Alanine[c] + peptidoglycan precursor[c] => D-Alanine[e] + peptidoglycan[c]'};
model=addRxns(model,rxns,3,'',true);                                                   % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'peptidoglycan'))),:)=[0];         % Assigns a confidence score of 0 as the reaction is artificial.


rxns.rxns={'peptidoglycan precursor biosynthesis'};                                                           % Creates this new reaction ID.
rxns.equations={'Undecaprenyl-diphospho-N-acetylmuramoyl-(N-acetylglucosamine)-L-alanyl-D-glutamyl-meso-2,6-diaminopimeloyl-D-alanyl-D-alanine => di-trans,poly-cis-Undecaprenyl diphosphate + peptidoglycan precursor'};   
% NB! It is unclear whether the metabolite with the long name in the reaction above is the correct compound!
model=addRxns(model,rxns,2,'c',true);                                                                         % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'peptidoglycan precursor biosynthesis'))),:)=[0];         % Assigns a confidence score of 0 as the reaction is artificial.


rxns.rxns={'protein'};                                                                  % Creates this new reaction ID.
rxns.equations={'0.512 L-Glutamine + 0.43 L-Phenylalanine + 0.687 L-Valine + 1.211 L-Alanine + 0.512 L-Glutamate + 1.135 Glycine + 0.115 L-Cysteine + 0.369 L-Aspartate + 0.159 L-Methionine + 0.369 L-Asparagine + 0.223 L-Histidine + 0.306 L-Isoleucine + 0.456 L-Arginine + 0.764 L-Threonine + 0.997 L-Proline + 0.189 L-Lysine + 40 ATP + 0.522 L-Leucine + 0.008 L-Tryptophan + 0.421 L-Serine + 0.222 L-Tyrosine => protein + 40 ADP'};
model=addRxns(model,rxns,2,'c',true);                                                   % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'protein'))),:)=[0];                % Assigns a confidence score of 0 as the reaction is artificial.


rxns.rxns={'cofactors_and_vitamines'};                                                  % Creates this new reaction ID.
rxns.equations={'0.167 NAD+ + 0.418 Thiamine + 0.14 Ubiquinone-8 + 0.141 FAD + 0.0074 Heme A + 0.249 Tetrahydrofolate + 0.243 FMN + 0.149 NADP+ + 0.145 CoA + 0.449 Pyridoxal phosphate => cofactors and vitamines'};
model=addRxns(model,rxns,2,'c',true);                                                   % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'cofactors_and_vitamines'))),:)=[0];% Assigns a confidence score of 0 as the reaction is artificial.


rxns.rxns={'DNA'};                                                                      % Creates this new reaction ID.
rxns.equations={'1.054 dCTP + 0.564 dTTP + 1.054 dGTP + 0.564 dATP + 4.4 ATP => 4.4 Orthophosphate + DNA + 4.4 ADP'};
model=addRxns(model,rxns,2,'c',true);                                                   % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'DNA'))),:)=[0];                    % Assigns a confidence score of 0 as the reaction is artificial.


rxns.rxns={'Biomass'};                                                                  % Creates this new reaction ID.
rxns.equations={'0.034 lipopolisaccharide + 0.06 RNA + 0.055 carbohydrate + 0.0495 phospholipid + 0.06 peptidoglycan + 0.68 protein + 0.03 cofactors and vitamines + 0.031 DNA + 15.3 ATP => Biomass + 15.3 Orthophosphate + 15.3 ADP'};
model=addRxns(model,rxns,2,'c',true);                                                   % Adds the new reaction to the model.    
model.rxnConfidenceScores((find(ismember(model.rxns,'Biomass'))),:)=[0];                % Assigns a confidence score of 0 as the reaction is artificial.


                                                                            % A transport and a producing exchange reaction is added for Biomass:

[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Biomass'},true,false);              % Adds a reversible transport reaction for Biomass.
model.rxns(end,1)={'TRP_c->e_Biomass'};                                                 % Assigns an appropriate reaction ID to the transport reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_Biomass'))),:)=[0];       % Assigns a confidence score of 0.
lastMet=model.mets(end);                                                    
[model, addedRxns]=addExchangeRxns(model,'out',lastMet);                                % Adds a producing exchange reaction for Biomass.
model.rxns(end,1)={'EXC_OUT_Biomass'};                                                  % Assigns an appropriate reaction ID to the exchange reaction.
model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_OUT_Biomass'))),:)=[0];        % Assigns a confidence score of 0.


                                                                            %% Validation/testing

model=setParam(model, 'obj', {'Biomass'},[1]);                              % The objective is set to maximize production of biomass.                                                                            
                                                                            
                                                                            % For the sake of error checking, constraints relating to energy and redox balance are
                                                                            % temporarily relaxed by the introduction of these fake reactions:
                                                                            
%rxns.rxns={'FAKE_FREE_ATP';'FAKE_FREE_NADH';'FAKE_FREE_NADPH'};             % Creates these new reaction IDs.
%rxns.equations={'ATP <=> ADP';'NAD+ <=> NADH';'NADP+ <=> NADPH'};           % Creates the new reaction equations.
%model=addRxns(model,rxns,2,'c');                                            % Adds the new reactions to the model.    
%model.rxnConfidenceScores((find(ismember(model.rxns,'FAKE_FREE_ATP'))),:)=[0];          % Assigns a confidence score of 0
%model.rxnConfidenceScores((find(ismember(model.rxns,'FAKE_FREE_NADH'))),:)=[0];         % Assigns a confidence score of 0
%model.rxnConfidenceScores((find(ismember(model.rxns,'FAKE_FREE_NADPH'))),:)=[0];        % Assigns a confidence score of 0
                                                    
                                                                            % In an additional preparation for error checking, fake transport and exchange reactions 
                                                                            % are temporarily added for H+ and Orthophosphate:
                                                                           
                                                                                                                                                        
%[model, addedRxns]=addTransport(model,{'c'},{'e'},{'H+'},true,false);                   % Adds a reversible transport reaction for H+.
%model.rxns(end,1)={'TRP_c->e_H+'};                                                      % Assigns an appropriate reaction ID to the transport reaction.
%model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_H+'))),:)=[0];            % Assigns a confidence score of 0.
%lastMet=model.mets(end);                                                    
%[model, addedRxns]=addExchangeRxns(model,'both',lastMet);                               % Adds a reversible exchange reaction for H+.
%model.rxns(end,1)={'EXC_OUT_H+'};                                                       % Assigns an appropriate reaction ID to the exchange reaction.
%model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_OUT_H+'))),:)=[0];             % Assigns a confidence score of 0.
                                                                            
%[model, addedRxns]=addTransport(model,{'c'},{'e'},{'Orthophosphate'},true,false);       % Adds a reversible transport reaction for Orthophosphate.
%model.rxns(end,1)={'TRP_c->e_Orthophosphate'};                                          % Assigns an appropriate reaction ID to the transport reaction.
%model.rxnConfidenceScores((find(ismember(model.rxns,'TRP_c->e_Orthophosphate'))),:)=[0];% Assigns a confidence score of 0.
%lastMet=model.mets(end);                                                    
%[model, addedRxns]=addExchangeRxns(model,'both',lastMet);                               % Adds a reversible exchange reaction for Orthophosphate.
%model.rxns(end,1)={'EXC_OUT_Orthophosphate'};                                           % Assigns an appropriate reaction ID to the exchange reaction.
%model.rxnConfidenceScores((find(ismember(model.rxns,'EXC_OUT_Orthophosphate'))),:)=[0]; % Assigns a confidence score of 0.



                                                                            % The fake reations for ATP, NADH and NADPH introduced earlier can now be removed:            
%model=removeReactions(model,{'FAKE_FREE_ATP';'FAKE_FREE_NADH';'FAKE_FREE_NADPH'});     % Removes the fake reactions introduced earlier.  

                                                                            %% Saving the GEM
    
exportToExcelFormat(model,'HPseGEM.xlsx')                                   % Saves the resulting GEM as an Excel-file (.xlsx-format) called HPseGEM.
exportModel(model,'HPseGEM.xml');                                           % Saves the resulting GEM as an SBML-file (.xml-format) called HPseGEM.
