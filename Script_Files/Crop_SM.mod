############################################
# Model:Estimation of Evapration( Simulated Model for Crop)                                #
# Date:15th July, 2017                                    #
# Version: 1.0                             #
# Author: Gifty Attiah 
##(References: Initial Data obtatined from Dr. Georg Hörmann(Digital Spatial Analysis)                                 #
############################################


binding
Raintimeseries=prec.tss; #Precipitation
pzones=zonee_rec.map;# zones for precipitation
ezones=zonee_rec.map;# zones for evaporation
ETptimeseries=etp.tss;  #Potential Evapration
##Landuse map and legend was edited from lu.map to Landuse.map to assign co-efficient values for diffrent landuse types
Landuse=landuse.map; #Landuse map with legend showing different landcover (jpeg file with legend)
outlet= dispoints_all.map; #Map showing all 8 oulets
#Discharge map created showing the highest simulated discharge and closest to soltfeld with map edits
gauging=discharge.map;
Soil=soil.map;# Map showing diffrent soils in Kielstau
slope = slope_rec.map; #Map showing the diffrent slope values
#kinematic function
 Beta=scalar(0.6);
 LDD=ldd_rec.map; # Map showing Local Drain Direction 
 N=Manning_rec.map;			# Manning's N;
 ChannelN=channelM_rec.map; # Map about channels
 ChannelWidth=channelw_rec.map; #Map showing width of Channels
 T=scalar(86400); # Time Scale for a day showing number of seconds

areamap
clone.map;

timer
1 9860 1 ;# Time scale output for 9860 days
reportdefault= 365 + 365 ..endtime;

initial
Inter_cap=scalar(1); #Interception capacity set to 1
soil_field_capacity =scalar(200.00); #Water content at field capacity
soil_start_reduction=scalar(100.00); # reduction point
Infil_fac=scalar(0.3); #Infiltration factor
SWC = scalar(0.001); #soil water content 
gw_factor = 0.1; # discharge flux from groundwater
gw_content = scalar(50); #Ground water content scale
#Proposed evaporation co-effcicient for Crop
Crop=scalar(1.2);

day = 0;

Q2_sum = scalar(0.001);

#kinematic function
 CL=celllength(); #Length of cell
 CA=cellarea(); #Area of cell
 DCL=max(downstreamdist(LDD),CL);
 ##channel width
 #? ChannelWidth=cover(ChannelWidth,0);
 ##initial streamflow(m3/s)
 Q=scalar(0.00000001);
 ##initial water height(m)
 H=scalar(0.000000001);

 # covers the 0 value in slope.map with 0.0001
 slope = if(slope eq 0, 0.0000001, slope);

 ##term for Alpha
 AlphaFact=(ChannelN)/(sqrt(slope))**Beta;
 #Power for Alpha; slope.map developed from dgm.map using slope operator
 AlphaPower=(2/3)*Beta;
 WH=0.000;
 FlowWidth= ChannelWidth;
 Q2=0.0001;

dynamic
Precipitation = timeinputscalar(Raintimeseries, pzones);# Measured values for precipitation
ETp = timeinputscalar(ETptimeseries, ezones)*Crop;# Potential evaporation 
#Estimated Evapotranspiration(ET)
EET_Crop.tss=timeoutput(pzones, ETp);

day = if (day >365 then 0 else day + 1);
#Calculating interception
interception_content = min(Precipitation, Inter_cap); 
interception = min(ETp, interception_content);#Inteception 

rest_et = ETp - interception;# rest of evaporation calculated by potential evaporation minus intercepted precipitation
rest_precipitation = Precipitation - interception;# Water to the soil before infiltration

Infiltration = (soil_field_capacity - SWC)*Infil_fac;	# infiltrated water
runoff_surface = max(rest_precipitation - Infiltration, 0);	# initial rapid surface runoff

rest_precipitation = rest_precipitation - runoff_surface; # remaining precipitation after rapid surface runoff

SWC=SWC+rest_precipitation; #remaining precipitation+already available water content
#Actual Evaporation for Crop
ETa = if (SWC > soil_start_reduction then rest_et else rest_et*(SWC/soil_start_reduction));	# actual evapotranspiration
ETa_Crop.tss=timeoutput(pzones,ETa);#actual evaporation for Crop output .tss file

#Soil water content for Crop
SWC = SWC - min(ETa, SWC);#soil water content after actual evaporation

Overflow = max(SWC-soil_field_capacity,0); #Initial overflow
gw_flux = if(SWC>soil_start_reduction then SWC*0.2 else 0);
gw=Overflow+gw_flux;# Ground water 

SWC=SWC-gw; #Final soil water content
SWC_Crop.tss=timeoutput(pzones,SWC); #output .tss for Soil water content for Crop

gw_content = gw_content + gw;
gw_outflow = gw_content * gw_factor; # Ground water overflow
gw_content = gw_content - gw_outflow;# Ground water content

# Runoff for Crop
Runoff = runoff_surface + gw_outflow;	# runoff is from both surface runoff and ground water outflow
Runoff_Crop.tss=timeoutput(gauging,Runoff);#output .tss for initial runoff for Crop

SumR=Runoff/1000;   #to m^3
   WH=WH+SumR;          # water height

 ## runoff calculation using one of kimematic
    ALPHA = AlphaFact*((FlowWidth+2*WH)**AlphaPower);
    QIn=SumR*CA/T/DCL;
#QIn and Q2 are both in m^3/S
Q2 = kinematic(LDD, Q2, QIn, ALPHA, Beta, 1,T,DCL);
    # error check from mailing list
    #Q2 = if(Q2<0.001,0.01,Q2);
    #Q2 = if(Q2>1e20,1,Q2);
    # Q2=cover(Q2,0.1);

    #flow velocity(m3/s)
    V=Q2/CA;   #CA is cell area
    #wh in mm unit
    WH=if(FlowWidth >0.001, 1000*(ALPHA*(Q2**Beta))/FlowWidth,0);

Q2_sum = Q2_sum + Q2;
# Output .tss kinematic runoff for all outlets for Crop
report Q2_Crop.tss=timeoutput(outlet,Q2);
## Output .tss kinematic runoff for virtual gauging station for Crop
report Q2_runoff_Crop.tss=timeoutput(gauging,Q2);