clear all
close all

%Parameters k and d are found from Max/Martin
%k s. Table S3.5 (Martin)
%d s. Table S3.6 (Martin)

%The birth parameter can be estimated from the number of filopodia births
%per growth cone and in one hour (at P40, values are left in comments) or
%analyticaly.

%c2 is a parameter for the tanh-function regulating the negative growth
%feedback over time. With the decay of born filopodia at P40 and P60 (in
%percent, x marks how much is still left at time t), evaluate
%c2 = t/(1+1/3*atanh(1-2*x)).

%Decay of stabilized filopodia is the inverse of their average life time.
%d_sF s. Table S3.09

%The bulbous stabilization parameters k_B are up till now equal for both
%filopodia and stabilized filopodia. That could be changed.
%The calculation behind it, follows from 
%r_->B*(B50/B50+B) = (r_F->B + r_sF->B)*(B50/B50+B)  = (k_F*F + k_sF*sF)*(B50/B50+B).
%setting both ks equal to each other (we'll name it k_F/sF) and rearranging, it follows that 
% k_F/sF = r_->B/(F+sF);
%Where we will substitute the average values at P60 for F and sF into the equation (Table S3.3 & S3.5). The parameter F still needs to be
%corrected with undetect. filo., i.e.: F(detected) = P(detected) * F,
%so in the end, we get:
% k_F/sF = r_->B/(F(detected)/P(detected)+sF)
%The reaction rate r_->B (r_2B) is fitted to the data (Table S.19)

%Decay rate of bulbous tips is again the inverse of the average life time.
%d_B s. Table S3.13

%The half maximum B_50 plays a role in the product-inhibition. If this
%parameter is high, the associated reaction will be activated very few times.

%The molecules array stores, which molecules are present (1) or deactivated
%(0).

%In the stabilization of filopodia we have a poisson shift due to the
%filopodiae length. The rate of stabilization (rateStab) and the hill
%parameter are calculated from its own routine.

%The parameters F_P40 and sF_P40 measured the average number of filopodia
%and stabilized filopodia at P40. The filopodia number needs to be
%corrected by the undetected ones.

param_set = {'T18C'; 'T25C';'T29C'};


for i=1:length(param_set)
    mutant = char(param_set(i));
    switch mutant
        
        case 'T18C'
            retract = 0;
            % Filopodia dynamics
            Parameters.(mutant).c1_sF = 1.823;%birth           
            Parameters.(mutant).c1_LF = 0.28;%birth
            Parameters.(mutant).c2_sF = 0.43;%death
            Parameters.(mutant).c2_LF = 0.07;%death
            %bulbous dynamics
            Parameters.(mutant).c3 = 0.3056;%0.0245;%;%;%rF2B
            Parameters.(mutant).B50 = 0.1187;%0.0522;%rB2S
            Parameters.(mutant).c4 = 0.0706;%1/30;%1/109;%rB20
            Parameters.(mutant).c5 = 0.0427;%0.0683;%;%rB2S
            Parameters.(mutant).c6 = 1/300;%1/150;%0.0683;%;%rB2S
            Parameters.(mutant).c7 = 0.006;%1
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.(mutant).p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = 1000;
            Parameters.(mutant).h = 1;
            
        case 'T25C'
             retract = 0;%0.85;
            % Filopodia dynamics
            Parameters.(mutant).c1_sF = 1.823;%birth           
            Parameters.(mutant).c1_LF = 0.28;%birth
            Parameters.(mutant).c2_sF = 0.43;%death
            Parameters.(mutant).c2_LF = 0.07;%death
            %bulbous dynamics
            Parameters.(mutant).c3 = 0.0860;%0.2956;%
            Parameters.(mutant).B50 = 0.3505;%0.0584;%
            Parameters.(mutant).c4 = 0.0706;%
            Parameters.(mutant).c5 = 0.0207;%0.2343;%
            Parameters.(mutant).c6 = 1/110;%1/100;%0.0683;%;%rB2S
            Parameters.(mutant).c7 = 0.006;%1
           % time dependent functions
            % a) dampening of filo dynamics
            Parameters.(mutant).p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = 1000;
            Parameters.(mutant).h = 1;

        case 'T29C'  
             retract = 0;%0.6;
            % Filopodia dynamics
            Parameters.(mutant).c1_sF = 1.823;%birth           
            Parameters.(mutant).c1_LF = 0.28;%birth
            Parameters.(mutant).c2_sF = 0.43;%death
            Parameters.(mutant).c2_LF = 0.07;%death
            %bulbous dynamics
            Parameters.(mutant).c3 = 0.0569;%0.1466;%0.0274;%
            Parameters.(mutant).B50 = 0.1044;%0.0489;%rB2S
            Parameters.(mutant).c4 = 0.0706;%0.185;% %rB20
            Parameters.(mutant).c5 = 0.0081;%0.3549;%rB2S
            Parameters.(mutant).c6 = 1/50;%1/60;%0.0683;%;%rB2S
            Parameters.(mutant).c7 = 0.006;%1
            % time dependent functions
            % a) dampening of filo dynamics
            Parameters.(mutant).p = [-2.9740e-17   3.3115e-13  -1.2896e-09   2.0637e-06  -1.4483e-03   1.0021e+00];
            % b) increase in propensity to form bulbous
            Parameters.(mutant).thalf = 1000;
            Parameters.(mutant).h = 1;

    end
end

try
    delete 'AllParametersSimple3.mat'
catch
end

save 'AllParametersSimple3.mat' Parameters 


