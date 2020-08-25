%Analytical floodplain overbank deposition model
%Written by: J. A. Nghiem
%Last edited: August 25, 2020

%Summary: This script runs the floodplain deposition model (see section
%Floodplain sedimentation model in supplement).

clear

%Model inputs, roughly scaled to Mississippi River
L=30000; %m, half-floodplain width to discretize spatial domain
l=1300; %distance over which to compute dimensionless floodplains statistics, take to be advection length of sand (m)
n=100000; %number of space steps to discretize floodplain
q=2.3; %m^2/s, overbank per-width discharge
sed_density=2650; %kg/m^3, sediment density
f_density=1000; %kg/m^3, fluid density (water)
g=9.81; %m/s^2, gravitational acceleration
t_nofloc=1; %time ratio to scale deposition rates to retrieve deposit thickness in un-flocculated case
t_floc=0.53; %time ratio to scale deposition rates to retrieve deposit thickness in flocculated case
%these time ratios ensure floodplain-integrated deposit thickness is
%identical in both flocculated and un-flocculated cases
ws_floc=[0.00006719 0.00034 0.001216]; %m/s, floc settling velocities
thresh_floc=[9.127E-06 2.053E-5 3.882E-5]; %m, flocculated grain size threshold (all sizes below are flocculated)
%ws_floc and thresh_floc are vectors of the same length representing
%different flocculation scenarios
%here, they represent lower, average, and upper scenarios for floc settling
%velocity based on estimates by Lamb et al., 2020

%Begin calculations below
%Read in grain size discretization and suspended sediment concentrations
%generated from the script parametric_gsd.R
sc=readtable('parametric_gsd.csv');
%Extract data from the table
d=sc{:,'center'}/1000; %particle diameter (m)
gsc=sc{:,'gsc'}; %grain size specific volumetric sediment concentration

x=linspace(0, L, n); % discretize floodplain into evenly spaced steps
R=(sed_density-f_density)/f_density; %submerged specific gravity of sediment

r0=ones(length(d), 1); %set a constant sediment concentration stratification of 1
lmb=zeros(1, length(x)); %set a constant deposit porosity of 0

%Initialize a multidimensional array to store model results
model_results_floc=NaN(length(d), length(x), length(ws_floc));

%Loop over each floc scenario and calculate grain size specific deposition
%rate as a function of distance from channel
for j=1:length(ws_floc)
    %Set the floc settling velocity
    wsfloc=ws_floc(j); %floc settling velocity (m/s)
    
    %Set the threshold for flocculated sizes
    bnd=thresh_floc(j); %convert from m to mm
    flocbnd=[0 bnd];
    
    %Calculate settling velocities as a function of grain size
    fsizes=(d>=flocbnd(1)) & (d<=flocbnd(2)); %find which grain sizes are flocculated
    ws=zeros(height(sc), 1); %initialize a vector to store settling velocities
    l_cent=d(~fsizes); %select un-flocculated grain sizes
    %Compute settling velocity using relation of Ferguson and Church, 2004
    %for un-flocculated sediment
    ws(~fsizes)=(R*g.*(l_cent.^2))./((20*1.0035*10^(-6))+sqrt(0.75*1.1*R*g*l_cent.^3));
    %Set uniform floc settling velocity
    ws(fsizes)=wsfloc;
    
    %Calculate products that appear in the model
    c1=ws.*r0.*gsc; %settling velocity, stratification, and sediment concentration
    c2=-ws.*r0/q; %settling velocity, stratification, divided by discharge
    
    %Initialize a matrix to store deposition rates
    %rows are different grain sizes
    %columns are different locations along the floodplain
    echange=NaN(length(d), length(x));
    
    for i=1:length(d) %loop over grain size classes
        echange(i,:)=c1(i).*exp(c2(i).*x);
    end
    %Apply porosity multiplicative factor
    echange=echange./(1-lmb);
    
    %Save results to cell array
    model_results_floc(:,:,j)=echange;
end

%Calculate results for un-flocculated case
l_cent=d; %select un-flocculated grain sizes
%Compute settling velocity using relation of Ferguson and Church, 2004
%for un-flocculated sediment
ws=(R*g.*(l_cent.^2))./((20*1.0035*10^(-6))+sqrt(0.75*1.1*R*g*l_cent.^3));

%Calculate products that appear in the model
c1=ws.*r0.*gsc; %settling velocity, stratification, and sediment concentration
c2=-ws.*r0/q; %settling velocity, stratification, divided by discharge

%Initialize a matrix to store deposition rates
%rows are different grain sizes
%columns are different locations along the floodplain
echange=NaN(length(d), length(x));

for i=1:length(d) %loop over grain size classes
    echange(i,:)=c1(i).*exp(c2(i).*x);
end
%Apply porosity multiplicative factor
echange=echange./(1-lmb);

%Save results to a matrix
model_results_nofloc=echange;

%Normalize deposition rates by maximum deposition rate out of all scenarios
%calculated
%combine floc and un-flocculated results
model_results_combined=cat(3, model_results_floc, model_results_nofloc);
%aggregate size classes by taking column sums
model_results_aggregated=sum(model_results_combined, 1);
%find the fastest deposition rate over all scenarios
max_rate=max(max(model_results_aggregated));
%normalize the results by the maximum rate
model_results_floc_normalized=model_results_floc./max_rate;
model_results_nofloc_normalized=model_results_nofloc./max_rate;
%rows are different grain sizes (in order from smallest to largest)
%columns are different locations (x=0 is at the channel, move away from
%channel as x increases)

%Multiply by time ratio factor to have same floodplain-integrated deposit
%thicknesses
model_results_floc_normalized=model_results_floc_normalized*t_floc;
model_results_nofloc_normalized=model_results_nofloc_normalized*t_nofloc;

%Calculate mud fraction (relative to sand)
mud_sizes=(d<=62.5*10^(-6)); %determine which grain sizes are mud (smaller than 62.5 microns)

%combine floc and un-flocculated results
model_results_combined=cat(3, model_results_floc_normalized, model_results_nofloc_normalized);
%aggregate size classes by taking column sums, represents total thickness
model_results_aggregated=sum(model_results_combined, 1);

%aggregate mud sizes
mud_aggregated=sum(model_results_combined(mud_sizes,:,:), 1);

%divide matrices to calculate mud fraction
mud_fraction_combined=mud_aggregated./model_results_aggregated;
%break out again into flocculated and un-flocculated cases
mud_fraction_floc=mud_fraction_combined(:,:,1:length(ws_floc));
mud_fraction_nofloc=mud_fraction_combined(:,:,end);

%Visualize mud fraction as a function of distance
figure
f1=plot(x, mud_fraction_floc(:,:,2), 'color', [126 47 142]./255, 'linewidth', 2); %flocculated
hold on
f2=plot(x, mud_fraction_nofloc, 'color', [0 114 189]./255, 'linewidth', 2); %un-flocculated
xlim([0 250])
ylim([0 1])
ylabel('mud fraction')
xlabel('distance from channel (m)')
legend([f1 f2], 'flocculated', 'un-flocculated', 'Location', 'northwest')