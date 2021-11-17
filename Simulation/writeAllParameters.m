clear all
close all

param_set = {'T18C'; 'T25C';'T29C'};


for i=1:length(param_set)
    mutant = char(param_set(i));
    switch mutant
        
        case 'T18C'
            retract = 0;
            % Filopodia dynamics
            Parameters.(mutant).null1 = 1;% empty reaction; for numerical accuracy         
            Parameters.(mutant).null2 = 0.5;% empty reaction; for numerical accuracy   
            %bulbous dynamics
            Parameters.(mutant).r1 = 0.3056;%0.0245;%;%;%rF2B
            Parameters.(mutant).B50 = 0.1187;%0.0522;%rB2S
            Parameters.(mutant).c2 = 0.0706;%1/30;%1/109;%rB20
            Parameters.(mutant).c3 = 0.0427;%0.0683;%;%rB2S
            Parameters.(mutant).c4 = 1/300;%1/150;%0.0683;%;%rB2S
            Parameters.(mutant).c5 = 0.006;%1
  
            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = 1000;
            Parameters.(mutant).h = 1;
            
        case 'T25C'
             retract = 0;%0.85;
            % Filopodia dynamics
            Parameters.(mutant).null1 = 1;% empty reaction; for numerical accuracy         
            Parameters.(mutant).null2 = 0.5;% empty reaction; for numerical accuracy   
            %bulbous dynamics
            Parameters.(mutant).r1 = 0.0860;%0.2956;%
            Parameters.(mutant).B50 = 0.3505;%0.0584;%
            Parameters.(mutant).c2 = 0.0706;%
            Parameters.(mutant).c3 = 0.0207;%0.2343;%
            Parameters.(mutant).c4 = 1/110;%1/100;%0.0683;%;%rB2S
            Parameters.(mutant).c5 = 0.006;%1

            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = 1000;
            Parameters.(mutant).h = 1;

        case 'T29C'  
             retract = 0;%0.6;
            % Filopodia dynamics
            Parameters.(mutant).null1 = 1;% empty reaction; for numerical accuracy         
            Parameters.(mutant).null2 = 0.5;% empty reaction; for numerical accuracy   
            %bulbous dynamics
            Parameters.(mutant).r1 = 0.0569;%0.1466;%0.0274;%
            Parameters.(mutant).B50 = 0.1044;%0.0489;%rB2S
            Parameters.(mutant).c2 = 0.0706;%0.185;% %rB20
            Parameters.(mutant).c3 = 0.0081;%0.3549;%rB2S
            Parameters.(mutant).c4 = 1/50;%1/60;%0.0683;%;%rB2S
            Parameters.(mutant).c5 = 0.006;%1
            
            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = 1000;
            Parameters.(mutant).h = 1;

    end
end

try
    delete 'AllParameters.mat'
catch
end

save 'AllParameters.mat' Parameters 


