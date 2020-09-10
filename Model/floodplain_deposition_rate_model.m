%Analytical floodplain overbank deposition rate model
%Written by: J. A. Nghiem
%Last edited: September 10, 2020

%Summary: This script runs the floodplain deposition model to calculate
%relative mud deposition rates (see Fig. 3D and the supplement for more
%details).

clear

%Inputs
max_L=10^6; %m, maximum distance from channel
n=1000; %number of log-spaced points to interpolate for distances from channel
q=2.3; %m^2/s, overbank per-width discharge
sed_density=2650; %kg/m^3, sediment density
g=9.81; %m/s^2, gravitational acceleration
por=0; %deposit porosity
%Set of floc settling velocities and their corresponding grain size
%flocculation thresholds
wsfloc=[0.00006719 0.00034 0.001216]; %m/s, floc settling velocities
bnd=[9.127E-06 2.053E-5 3.882E-5]; %m, flocculated grain size threshold (all sizes below are flocculated)
%wsfloc and bnd are vectors of the same length representing
%different flocculation scenarios
%here, they represent lower (16th percent quantile), average, and upper (84th percent quantile) scenarios for floc settling
%velocity based on estimates by Lamb et al., 2020
cutoff=62.5*10^(-6); %m, calculate deposition rates for grain sizes below this cutoff (mud)


%Begin calculations below
%Log interpolate floodplain lengths from 1 m to the specified maximum
L=exp(linspace(log(1), log(max_L), n));

%Read in the parametric grain size and suspended sediment concentration data
sc=readtable('parametric_gsd.csv');

%Extract data from the table
d=sc{:,'center'}; %mm, particle diameter
gsc=sc{:,'gsc'}; %grain size specific volumetric sediment concentration

%Initialize vectors to store floodplain deposition rate as function of
%distance from channel
dr_floc=NaN(1, length(L)); %flocculated case at average floc settling velocity
dr_floc_lower=NaN(1, length(L)); %flocculated case at lower floc settling velocity
dr_floc_upper=NaN(1, length(L)); %flocculated case at upper floc settling velocity
dr_nofloc=NaN(1, length(L)); %un-flocculated case

%Calculate settling velocities for each case
cutoff=cutoff*10^3; %convert cutoff from m to mm 
bnd=bnd*1000; %convert bnd from m to mm
cut_crit=(d<cutoff); %define cutoff criterion
d=d(cut_crit); %keep only grain sizes smaller than cutoff
gsc=gsc(cut_crit); %keep corresponding sediment concentration
R=(sed_density-1000)/1000; %submerged specific gravity of sediment
l_cent=d/1000; %select non-floc sizes and convert to m

%Compute settling velocity using Ferguson and Church (2004) for
%un-flocculated sediment
ws=(R*g.*(l_cent.^2))./((20*1.0035*10^(-6))+sqrt(0.75*1.1*R*g*l_cent.^3));
ws_flocculated=ws;
ws_flocculated_lower=ws;
ws_flocculated_upper=ws;

%Set uniform settling velocity for floc for lower, middle, and upper
%for grain sizes below each corresponding flocculation threshold
ws_flocculated(d<=bnd(2))=wsfloc(2);
ws_flocculated_lower(d<=bnd(1))=wsfloc(1);
ws_flocculated_upper(d<=bnd(3))=wsfloc(3);

r0=ones(length(d), 1); %set a constant sediment concentration stratification of 1

%Loop over different distances on the floodplain and calculate deposition
%rate for each flocculation case
for j=1:length(L)
    L_fp=L(j); %extract distance from the channel
    %Calculate mud deposition rate as a function of distance from channel
    %for different flocculation scenarios
    dr_floc(j)=sum((ws_flocculated.*r0.*gsc).*exp(-ws_flocculated.*r0*L_fp/q))./(1-por);
    dr_floc_lower(j)=sum((ws_flocculated_lower.*r0.*gsc).*exp(-ws_flocculated_lower.*r0*L_fp/q))./(1-por);
    dr_floc_upper(j)=sum((ws_flocculated_upper.*r0.*gsc).*exp(-ws_flocculated_upper.*r0*L_fp/q))./(1-por);
    dr_nofloc(j)=sum((ws.*r0.*gsc).*exp(-ws.*r0*L_fp/q))./(1-por);
end

%Normalize deposition rates by the maximum rate in the un-flocculated case
maxv=max(dr_nofloc); %find the maximum rate in the un-flocculated case
%Normalize all the rates by this maximum value
dr_floc=dr_floc/maxv;
dr_floc_lower=dr_floc_lower/maxv;
dr_floc_upper=dr_floc_upper/maxv;
dr_nofloc=dr_nofloc/maxv;

%Plot the results
%Plot colors
floc_color=[126 47 142]./255; %color representing flocculated case
nofloc_color=[0 114 189]./255; %color representing un-flocculated case

figure
area(L, dr_floc_upper, 'FaceColor', floc_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on
area(L, dr_floc_lower, 'FaceColor', 'white', 'EdgeColor', 'none');
f1=plot(L, dr_floc, 'color', floc_color, 'linewidth', 2);
f2=plot(L, dr_nofloc, 'color', nofloc_color, 'linewidth', 2);
set (gca, 'Xscale', 'log');
xlabel('distance from channel (m)')
ylabel('normalized mud deposition rate')
legend([f1 f2], 'flocculated', 'un-flocculated')
%shaded area represents ranges of possible model results within the ranges
%of plausible floc settling velocities