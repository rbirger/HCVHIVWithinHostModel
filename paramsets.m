%% Set up parameter sets

%Parameters that don't change
params.sh = 9;
params.alph = 5e-3;%ah_orig;1.75e-3;%
params.gamm = 2e-8;%ac_orig;1e-9;%
params.beta_h = 4.1e-7;
params.dh = 9e-3;
params.VL = 0;
 

%bistable param set - for SVR

params.sc=4365; %recruitment rate of Healthy Hepatocytes
params.r1=2.7; %density-dependent proliferation rate for healthy hepatocytes .31
params.c= 10.06; %clearance rate of HCV virus orig 11.5
params.p = 13.48; %rate of viral production per infected hepatocyte 3.48
params.Tcmax =4.016e6; %Hepatocyte carrying capacity
params.dc=1.06e-3; %death rate of healthy hepatocytes
params.beta_c = 7.3e-8; %infection rate of healthy hepatocytes orig 7.4
params.r2 = 7.52; %proliferation rate of infected hepatocytes orig 4.4 7.52 for low ac/ah4.52 9.52
params.delta_c_orig = .72 ;%-ah*sh/dh; %death rate of infected hepatocytes orig 2.765 .17 .62 .305 (for low ac/ah) .82

paramsALL = params;

ind_bistable = 1;

% SVR parameter set that is not bistable

params.sc=9109; %recruitment rate of Healthy Hepatocytes
params.r1=2.6; %density-dependent proliferation rate for healthy hepatocytes .31 2.8759
params.c= 7.27; %clearance rate of HCV virus orig 11.5 19.5061
params.p = 1.48; %rate of viral production per infected hepatocyte 12.1854 41.1
params.Tcmax =11.2e6; %Hepatocyte carrying capacity 9598853
params.dc=1.11e-3; %death rate of healthy hepatocytes 1.0104e-3
params.beta_c = 6.09e-7; %infection rate of healthy hepatocytes orig 9.5902e-7
params.r2 = 4.35; %proliferation rate of infected hepatocytes orig 4.0586
params.delta_c_orig = .73;%-ah*sh/dh; %death rate of infected hepatocytes 2.5702 1.9 2.44


paramsALL(2) = params;
ind_nonbist = 2;

%Null responder param set
params.sc=11.6e3; %recruitment rate of Healthy Hepatocytes
params.r1=.34; %density-dependent proliferation rate for healthy hepatocytes .31
params.c=16.9; %clearance rate of HCV virus orig 11.5
params.p = 10; %rate of viral production per infected hepatocyte
params.Tcmax = 9.6e6; %Hepatocyte carrying capacity
params.dc=1.37e-3; %death rate of healthy hepatocytes
params.beta_c = 9.6e-7; %infection rate of healthy hepatocytes orig 7.4
params.r2 = 1.8; %proliferation rate of infected hepatocytes orig 4.4
params.delta_c_orig = 1.806;%-ah*sh/dh; %death rate of infected hepatocytes orig 2.


paramsALL(3) = params;
ind_NullResp = 3;

%Partial Responder
params.sc=903; %recruitment rate of Healthy Hepatocytes
params.r1=1.7; %density-dependent proliferation rate for healthy hepatocytes .31 2.8759
params.c=16.1; %clearance rate of HCV virus orig 11.5 19.5061
params.p = 18.09; %rate of viral production per infected hepatocyte 12.1854
params.Tcmax =11.2e6; %Hepatocyte carrying capacity 9598853
params.dc=1.25e-3; %death rate of healthy hepatocytes 1.0104e-3
params.beta_c = 4.6e-7; %infection rate of healthy hepatocytes orig 9.5902e-7
params.r2 = 1.77; %proliferation rate of infected hepatocytes orig 4.0586
params.delta_c_orig = 2.009;%-ah*sh/dh; %death rate of infected hepatocytes 2.5702

paramsALL(4) = params;

ind_Part = 4;

%Relapse
params.sc=8096; %recruitment rate of Healthy Hepatocytes
params.r1=2.88; %density-dependent proliferation rate for healthy hepatocytes .31 2.8759
params.c=19.52; %clearance rate of HCV virus orig 11.5 19.5061
params.p = 12.2; %rate of viral production per infected hepatocyte 12.1854
params.Tcmax =9.6e6; %Hepatocyte carrying capacity 9598853
params.dc=1.01e-3; %death rate of healthy hepatocytes 1.0104e-3
params.beta_c = 9.6e-7; %infection rate of healthy hepatocytes orig 9.5902e-7
params.r2 = 4.06; %proliferation rate of infected hepatocytes orig 4.0586
params.delta_c_orig = 2.10;%-ah*sh/dh; %death rate of infected hepatocytes 2.5702

paramsALL(5) = params;

ind_Relapse = 5;

%Breakthrough

%HCV params - %Breakthrough
params.sc=8096; %recruitment rate of Healthy Hepatocytes
params.r1=2.88; %density-dependent proliferation rate for healthy hepatocytes .31 2.8759
params.c=19.52; %clearance rate of HCV virus orig 11.5 19.5061
params.p = 12.2; %rate of viral production per infected hepatocyte 12.1854
params.Tcmax =9.6e6; %Hepatocyte carrying capacity 9598853
params.dc=1.01e-3; %death rate of healthy hepatocytes 1.0104e-3
params.beta_c = 9.6e-7; %infection rate of healthy hepatocytes orig 9.5902e-7
params.r2 = 4.06; %proliferation rate of infected hepatocytes orig 4.0586
params.delta_c_orig = 2.10;%-ah*sh/dh; %death rate of infected hepatocytes 2.5702

paramsALL(6) = params;
ind_Breakthrough = 6;

%Neumann 1 
params.sc=14351; %recruitment rate of Healthy Hepatocytes
params.r1=2.09; %density-dependent proliferation rate for healthy hepatocytes .31
params.c= 23.5; %clearance rate of HCV virus orig 11.5
params.p = 20.03; %rate of viral production per infected hepatocyte
params.Tcmax =12.4e6; %Hepatocyte carrying capacity
params.dc=1.11e-3; %death rate of healthy hepatocytes
params.beta_c = 9.02e-7; %infection rate of healthy hepatocytes orig 7.4
params.r2 = 2.09; %proliferation rate of infected hepatocytes orig 4.4
params.delta_c_orig = 2.28 ;%-ah*sh/dh; %death rate of infected hepatocytes orig 2.765

paramsALL(7) = params;
ind_Neumann1 = 7;

%Neumann 2
params.sc=6377; %recruitment rate of Healthy Hepatocytes
params.r1=2.07; %density-dependent proliferation rate for healthy hepatocytes .31
params.c= 9.75; %clearance rate of HCV virus orig 11.5
params.p = 33.2; %rate of viral production per infected hepatocyte
params.Tcmax =4.2e6; %Hepatocyte carrying capacity
params.dc=1.37e-3; %death rate of healthy hepatocytes
params.beta_c = 1.5e-7; %infection rate of healthy hepatocytes orig 7.4
params.r2 = 1.3; %proliferation rate of infected hepatocytes orig 4.4
params.delta_c_orig = 1.78 ;%-ah*sh/dh; %death rate of infected hepatocytes orig 2.765

paramsALL(8) = params;
ind_Neumann2 = 8;
