%Analytical floodplain-averaged overbank deposition rate model
%Written by: J. A. Nghiem
%Last edited: August 25, 2020

%Summary: This script runs grain size specific floodplain-averaged
%deposition rate model for flocculated and un-flocculated cases. At each
%floodplain width and grain size, the script calculates deposition rate
%averaged over the entire floodplain.

clear

%Inputs
min_L=1/2; %m, minimum half-floodplain width
max_L=(10^6)/2; %m, maximum half-floodplain width
n=100; %number of points to interpolate between the minimum and maximum
q=2.3; %m^2/s, overbank discharge
sed_density=2650; %kg/m^3, sediment density
f_density=1000; %kg/m^3, fluid density (water)
g=9.81; %m/s^2, gravitational acceleration
por=0; %set porosity as a constant 0
ho=3; %m, overbank water depth
ws_floc=[0.00006719 0.00034 0.001216]; %m/s, range of floc settling velocities (lower, middle, upper)
thresh_floc=[9.127E-06 2.053E-5 3.882E-5]; %m, flocculated grain size threshold (all sizes below are flocculated)cutoff=62.5*10^(-3); %set cutoff in mm, ignore grain size above this cutoff
norm=boolean(1); %whether you want to normalize the deposition rate or not by maximum

%Calculations below
L=exp(linspace(log(min_L), log(max_L), n)); %interpolate points along floodplain

%Read in grain size and suspended sediment concentration data
sc=readtable('parametric_gsd.csv');

%Extract data from the table
d=sc{:,'center'}/1000; %m, particle diameter
gsc=sc{:,'gsc'}; %grain size specific volumetric sediment concentration

%Initialize vectors to store floodplain-averaged deposition rate
%functions of floodplain width
dr_floc=NaN(1, length(L)); %average floc settling velocity case
dr_floc_lower=NaN(1, length(L)); %low floc settling velocity case
dr_floc_upper=NaN(1, length(L)); %high floc settling velocity case
dr_nofloc=NaN(1, length(L)); %un-flocculated case

%Calculate grain size specific settling velocities
cut_crit=(d<62.5*10^(-6)); %define mud cutoff criterion
d=d(cut_crit); %keep only grain sizes smaller than cutoff
gsc=gsc(cut_crit); %keep corresponding sediment concentration
R=(sed_density-f_density)/f_density; %submerged specific gravity of sediment
%Compute settling velocity using relation of Ferguson and Church, 2004 for non-floc
ws=(R*g.*(d.^2))./((20*1.0035*10^(-6))+sqrt(0.75*1.1*R*g*d.^3));
ws_flocculated=ws;
ws_flocculated_lower=ws;
ws_flocculated_upper=ws;
%Set uniform floc settling velocities for lower, average, and upper floc
%scenarios
ws_flocculated(d<=thresh_floc(2))=ws_floc(2);
ws_flocculated_lower(d<=thresh_floc(1))=ws_floc(1);
ws_flocculated_upper(d<=thresh_floc(3))=ws_floc(3);

r0=ones(length(d), 1); %set a constant sediment concentration stratification of 1

%Calculate floodplain-averaged mud deposition rate for each floodplain
%width
for j=1:length(L)
    L_fp=L(j); %extract the floodplain width
    %Calculate the floodplain-averaged deposition rate
    dr_floc(j)=(q/(L_fp*(1-por)))*sum(gsc.*(1-exp(-(ws_flocculated.*r0)*L_fp/q)));
    dr_floc_lower(j)=(q/(L_fp*(1-por)))*sum(gsc.*(1-exp(-(ws_flocculated_lower.*r0)*L_fp/q)));
    dr_floc_upper(j)=(q/(L_fp*(1-por)))*sum(gsc.*(1-exp(-(ws_flocculated_upper.*r0)*L_fp/q)));
    dr_nofloc(j)=(q/(L_fp*(1-por)))*sum(gsc.*(1-exp(-ws.*r0*L_fp/q)));
end

%Plot the results of mud deposition rate as a function of floodplain width
L_width=L*2; %multiply by 2 to obtain full floodplain width

%Normalize deposition rates by the maximum deposition rate if specified
if norm
    %Calculate the maximum rate
    maxv=max([dr_floc dr_nofloc dr_floc_upper dr_floc_lower]);
    %Normalize rates for each case
    dr_floc=dr_floc/maxv;
    dr_floc_lower=dr_floc_lower/maxv;
    dr_floc_upper=dr_floc_upper/maxv;
    dr_nofloc=dr_nofloc/maxv;
    
    %vertical axis label
    vlab='Relative floodplain-averaged mud deposition rate';
else %or don't normalize
    vlab='Floodplain-averaged mud deposition rate (m/s)';
end

figure
area(L_width, dr_floc_upper, 'FaceColor', [126 47 142]./255, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on
area(L_width ,dr_floc_lower, 'FaceColor', 'white', 'EdgeColor', 'none');
f1=plot(L_width, dr_floc, 'color', [126 47 142]./255, 'linewidth', 2);
f2=plot(L_width, dr_nofloc, 'color', [0 114 189]./255, 'linewidth', 2);
set (gca, 'Xscale', 'log');
xlabel('Floodplain width (m)')
ylabel(vlab)
legend([f1 f2], 'flocculated', 'un-flocculated')