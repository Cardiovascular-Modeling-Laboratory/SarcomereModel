%{
const_gen.m

Created: 11/15/2020
Last Updated: 05/09/2023

written by Andrew Schmidt
Cardiomyocyte Laboratory
Henry Samueli School of Engineering
University of California, Irvine, Irvine, CA 92697

The purpose of this code is to generate the constants used in the SarcContractionCycle

        
Output: .mat file of constants for the half-sarcomere contraction simulation
%}

function const_gen()

%% Define constants

%Define spring constants
F = 0; %external force (-) at z-line and (+) at m-line
a_z = 0; %z-line node location when az is fixed
km = 6060e-3; %6060e-3; %spring constant for myosin filament; 6060 pN/nm in N/m
ka = 5230e-3; %5230e-3; %spring constant for actin filament; 5230 pN/nm in N/m
kxb = 3e-3; %3e-3; %spring constant for XB; 3e-3 (Tanner, 2012) (range: 1-5, Duke 1.3) pN/nm in N/m
k_titin = 10e-3; %spring constant for titin; 10 pN/nm in N/m
kF = 0.5e-3; %0.5e-3; %spring constant of external force spring (representing PDMS); 50 pN/nm in N/m************

%Define rate constants' constants for state transitions 
T = 300; %in K
kB = 1.380649e-23; %in J/K
B = 1/(kB*T); %1/kBT in J^-1 (=2.4143e+20)
dG_bind = -4*kB*T; %in J
dG_stroke = -4.5*kB*T; %in J
dG_hyd_st = 13*kB*T; %in J
dG_c = -14*kB*T; %in J, the energy associated with dissociation of ADP and Pi during power stroke

% Define normalization concentrations
ATP_norm = 1; %in M
ADP_norm = 1; %in M
P_norm = 1; %in M
dG_hyd = kB*T*log((ATP_norm)/((P_norm)*(ADP_norm))) + dG_hyd_st; %in J
norm_by = '1M'; 
conc_factor = 1;
conc_factor_ADP = conc_factor; 
conc_factor_ATP = conc_factor; 
conc_ADP_actual = conc_factor_ADP*0.03e-3; %in M
conc_P_actual = 3e-3; %in M
conc_ATP_actual = conc_factor_ATP*5e-3; %in M
conc_ADP = conc_ADP_actual/ADP_norm; 
conc_P = conc_P_actual/P_norm; 
conc_ATP = conc_ATP_actual/ATP_norm;

kbind = 800; %in s^-1, determines rate of binding when there is no cross bridge deformation i.e. when actin and myosin xb head are perfectly aligned
k23_cap = 300; %in s^-1, rate of power stroke when there is no cross bridge deformation [state 2->3]
kATP0 = 1e-2; %in s^-1, rate of reverse hydrolysis and ATP dissociation from myosin cross bridge [state 1->3]

%Define constants for calcium regulation of actin availability
conc_Ca = 0; %initialize as zero calcium in system; 0 = Ca present, 1 = Ca not present
t_Ca_s = 0; % time at which calcium is introduced into the system
t_Ca = 1; % iteration at which calcium is introduced into the system
t_Ca_counter = 0; %initialize counter outside of "for loop"
Ca_reg = 0; %To ensure input of Ca "if-statement" only occurs once inside "for loop"

%% Save Constants
%%{           
filename = ['const_gen_',...
            'km=', num2str(km), ',', ...
            'ka=', num2str(ka), ',', ...
            'kxb=', num2str(kxb), ',', ...
            'k_titin=', num2str(k_titin), ',', ...
            'kF=', num2str(kF), ',', ...
            'conc_ADP=', num2str(conc_ADP), ',', ...
            'conc_P=', num2str(conc_P), ',', ...
            'conc_ATP=', num2str(conc_ATP), ',', ...
            'norm_by=', norm_by, ',', ...
            'dG_bind=', num2str(dG_bind*B), ',', ...
            'dG_stroke=', num2str(dG_stroke*B), ',', ...
            'dG_c=', num2str(dG_c*B), ',', ...
            'kbind=', num2str(kbind), ',', ...
            'k23cap=', num2str(k23_cap), ',', ...
            'kATP0=', num2str(kATP0), '.mat'];
%Template to add more lines              
%'=', num2str(), ',', ...
                    
save(filename);
%}








