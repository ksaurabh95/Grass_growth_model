# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 10:45:28 2025

@author: saurabh kumar
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 12:32:16 2025

@author: saurabh kumar
"""

import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(precision=17)
# soil type loam soil 
om= 1.724*2.4
sand = 44
silt = 32
clay = 24
SatCondDay = 98.2
InfRate = 360   # notes 
depth = 100
Nmin = 0.29
Norg =  0.29
To_mean = 18
maxT = 25


# Constants and initial values
soil_type_data = {'SatCondDay': SatCondDay, 'InfRate': InfRate}  # sandy soil
constants = {'propRunOff': 0.8, 'Ip': 0.004, 'K': 0.115, 'TrefMin': 15, 'TrefIm': 10, 'c': 0.2}
assumed_values = { 'depth': depth, 'OM': om, 'clay': clay, 'sand': sand, 'Norg': Norg, 'Nmin': Nmin,
                  'To_mean': To_mean, 'maxT': maxT   }

# Weather assumptions for 30 days


# importing meterological data
df = pd.read_excel('weather_data_2000_2025.xlsx')  # For .xlsx files

df['date'] = pd.to_datetime(df['date'])

rainfall = np.array(df["rain"])  # mm/day
ETP = np.array(df["pe"])   # mm/day


# Soil property calculations
def calc_water_capacity(sand, clay, OM, depth):
    waterCapacity = (0.2576 - 0.002 *sand + 0.0036 * clay + 0.0299 * OM) * pow(10,7) * depth / 100
    return waterCapacity 

def calc_saturation(waterCapacity):
    waterSaturation = 100 / 88 * waterCapacity
    return waterSaturation

def calc_wilting_point(clay, OM, depth):
    wiltingPoint = (0.026 + 0.005 * clay + 0.0158 * OM) * pow(10,7) * depth / 100
    return wiltingPoint






#--------------- Water model functions ------------------
    
def calc_water1( water_prev, rain, ETP, notRunOff_prev, InfRate):
    # ETP = evapotranspiration(mm/day) 
    if rain <= InfRate:
        water1 = water_prev + (rain - ETP )*10000 + notRunOff_prev
    else: 
        water1 = water_prev + (InfRate - ETP )*10000 + notRunOff_prev
        
    return water1

  
def calc_water_extra(rain,InfRate):
    # Depending on the infiltrationrate relative to the amount of rain received 
    # eamount of water that cannot be drained
    if InfRate < rain:
        waterExtra = rain - InfRate
    else: 
        waterExtra = 0
    
    return waterExtra


def calc_not_runoff(water1, Saturation, waterExtra, propRunOff):
    # Depending on the infiltrationrate relative to the amount of rain received 
    # will either runoff or stay at the paddock surface in mm
    if water1 > Saturation :
        notRunOff = (1-propRunOff)*(water1 - Saturation + waterExtra )
    elif waterExtra > 0:
        notRunOff = waterExtra*10000*(1-propRunOff)
    else:
        notRunOff = 0
    
    return notRunOff 



def calc_NoneSatConDay(water1,capacity, saturation, SatCondDay):
    if capacity <= water1 <= saturation:
        NoneSatConDay = SatCondDay*np.exp(48.2 * (water1 - saturation)*pow(10,-7) )
    else:
        NoneSatConDay = 0
    
    return NoneSatConDay
        

def calc_drained(water1, capacity, saturation, SatCondDay):
    if water1 < capacity:
        drained = 0
    elif capacity <= water1 <= saturation:
        NoneSatConDay = calc_NoneSatConDay(water1, capacity, saturation, SatCondDay)
        drained = min(NoneSatConDay * (water1 - capacity) / 100, water1 - capacity)
    else:
        drained = min(SatCondDay * (saturation - capacity) / 100, saturation - capacity)
    
    return drained


def calc_water(water1, drained, capacity, saturation):
    if water1 > capacity:
        water = water1 - drained
    elif water1 > saturation:
        water = saturation - drained
    else: water =  water1
    
    return water 


def calc_water10(water10_prev, rain, ETP, notRunOff_prev, wilting, saturation, depth):
    return max((wilting * 10 / depth),
               min(water10_prev + (rain - ETP) * 10000 + notRunOff_prev,
                   (saturation * 10 / depth)))

def calc_water_stress_depth(water, capacity, wilting):
    if water < wilting:
        return 0
    return (water - wilting) / (capacity - wilting)

def calc_water_stress_10cm(water10, wilting, capacity, depth):
    return (water10 - wilting * 10 / depth) / ((capacity - wilting) * 10 / depth)

def calc_Fw_phi(wstress, maxT):
    return (-1.2387 * wstress**2 + 2.2387 * wstress - 0.0056) * 18 / maxT


#------------------ mineralisation section -------------------------------


def calc_organic_Nitrogen(OM):
    if OM >3:
        N_org_corr = (3+ 6 * (1 - np.exp(-0.22*(OM-3)) ))*4200/1.75
    else: N_org_corr = OM*4200/1.75
    return N_org_corr

def calc_potential_mineralisation_rate(N_org_corr):
    Vp_min = 0
    Vp_max = 2.2
    Vp = ( 0.0929 + (0.1833-0.0929)*np.exp(-0.2173*N_org_corr/1000) )*N_org_corr/1000
    Vp = np.clip(Vp, Vp_min, Vp_max)
    
    
    return Vp


def calc_Temp_para_mineralisation(K, To_mean, TrefMin):
    if To_mean<0:
        fT_m = 0
    elif 0<=To_mean < 4:
        fT_m = To_mean*0.28/4 
    else:fT_m = np.exp(K*(To_mean - TrefMin))
    
    return fT_m



def calc_water_para_mineralisation(c, water, capacity, wilting):
    g0_min = c
    g0_max = 2.2
    g0 = (1-c)*water*wilting/(capacity*wilting) + c
    g0 = np.clip(g0, g0_min, g0_max)
    
    
    return g0


def calc_mineralisation(K, To_mean, TrefMin, c, water, capacity, wilting, OM ):
    g0 = calc_water_para_mineralisation(c, water, capacity, wilting)
    fT_m = calc_Temp_para_mineralisation(K, To_mean, TrefMin)
    N_org_corr = calc_organic_Nitrogen(OM) 
    Vp = calc_potential_mineralisation_rate(N_org_corr)
    
    
    return Vp*fT_m*g0

# -------------------------------- Immobilisation rate ------------------------




def calc_Temp_para_immobilisation(K, To_mean, TrefMin):
    if To_mean<0:
        fT1 = 0
    elif 0<=To_mean < 4:
        fT1 = To_mean*2/4 
    else:fT1 = np.exp(-K*(To_mean - TrefMin))
    
    return fT1


def calc_immobilisation(Ip, Nmin , c, water, capacity, wilting, OM ):
    g0 = calc_water_para_mineralisation(c, water, capacity, wilting)
    if Nmin <= 60:
        Immobilisation_rate = Ip*Nmin*g0*2.5 
    else: 
        Immobilisation_rate = Ip*60*Nmin*g0*2.5 
    
    return Immobilisation_rate


# -------------------------------- soil Nitrification and denitrification ------------------------


def calc_repartition_N2_N20(clay):
    repartition_N2_N20 = 0.189 + (1.71*clay)*0.01/ (1 + 1.36*clay*0.01) 
    return repartition_N2_N20



def calc_emission_N2_N20( Nmin,g0,repartition_N2_N20,fT_m ):
    
    repartition_N2_N20 = calc_repartition_N2_N20(clay)
    globalemmiosn = Nmin*fT_m*g0/1000
    N2 = (1-repartition_N2_N20)*globalemmiosn
    N20 = repartition_N2_N20*globalemmiosn
    
    return N2,N20








# Initial values
depth = assumed_values['depth']
OM = assumed_values['OM']
clay = assumed_values['clay']
sand = assumed_values['sand']

capacity = calc_water_capacity(sand, clay, OM, depth)
saturation = calc_saturation(capacity)
wilting = calc_wilting_point(clay, OM, depth)

water_prev = saturation
water10_prev = saturation*10/depth
notRunOff_prev = 0

results = []

for day in range(len(rainfall)):
# for day in range(100):

    rain = rainfall[day]
    etp = ETP[day]
    water1 = calc_water1(water_prev, rain, etp, notRunOff_prev, soil_type_data['InfRate'])
    water_extra = calc_water_extra(rain, soil_type_data['InfRate'])
    not_runoff = calc_not_runoff(water1, saturation, water_extra, constants['propRunOff'])
    drained = calc_drained(water1, capacity, saturation, soil_type_data['SatCondDay'])
    water = calc_water(water1, drained, capacity, saturation)
    water10 = calc_water10(water10_prev, rain, etp, notRunOff_prev, wilting, saturation, depth)
    stress_depth = calc_water_stress_depth(water, capacity, wilting)
    stress_10cm = calc_water_stress_10cm(water10, wilting, capacity, depth)
    stress = max(stress_depth, stress_10cm)
    Fw_phi = calc_Fw_phi(stress, assumed_values['maxT'])
    
    N_org_corr = calc_organic_Nitrogen(OM)
    Vp = calc_potential_mineralisation_rate(N_org_corr)
    fT_m = calc_Temp_para_mineralisation(constants['K'], assumed_values['To_mean'],constants['TrefMin'])
    g0 = calc_water_para_mineralisation(constants['c'], water, capacity, wilting)
    mineralisation_rate = calc_mineralisation(constants['K'], assumed_values['To_mean'], constants['TrefMin'], constants['c'], water, capacity, wilting, assumed_values['OM'] )
    
    fT1 = calc_Temp_para_immobilisation(constants['K'], assumed_values['To_mean'], constants['TrefMin'])
    Immobilisation_rate = calc_immobilisation(constants['Ip'], assumed_values['Nmin'] , constants['c'], water, capacity, wilting,assumed_values['OM'])

    repartition_N2_N20 = calc_repartition_N2_N20(assumed_values['clay']) 
    
    
    N2,N20 = calc_emission_N2_N20( Nmin,g0,repartition_N2_N20,fT_m )
    
    
    results.append([day+1,rain,etp, water1, water_extra, not_runoff, drained, water, water10, stress_depth, stress_10cm, stress, Fw_phi,
                    N_org_corr, mineralisation_rate, Immobilisation_rate,repartition_N2_N20,N2,N20])
    
    water_prev = water
    water10_prev = water10
    notRunOff_prev = not_runoff

# Save to file
sim_df = pd.DataFrame(results, columns=['Day','rain','etp', 'water1', 'waterExtra', 'notRunOff', 'drained', 'water', 
    'water10', 'waterStress_depth', 'waterStress_10cm', 'waterStress_used', 'Fw_phi','N_org_corr', 'mineralisation_rate',
    'Immobilisation_rate','repartition_N2_N20','N2','N20'])
sim_df["date"] = df["date"]
print(sim_df)
start_date = "2010-01-01"
end_date = "2020-01-01"  # Adjust if you get data beyond 1995

plt.figure(figsize=(12,4))
# plt.plot(sim_df["date"],sim_df["water10"]/10000,label = 'top 10 cm')
plt.plot(sim_df["date"],sim_df["water"]/10000,label = '100 cm')
plt.plot(sim_df["date"],sim_df["rain"],label = 'rain')
plt.plot(sim_df["date"],sim_df["etp"],label = 'etp')

plt.ylabel("Water")
# plt.xlabel("Day")
# plt.xlim([0,2000])
plt.legend()
# plt.yticks(range(400, 1201, 200))
plt.grid(True)
plt.xticks(pd.date_range(start=start_date, end=end_date, freq="YS"), rotation=45)  # 'YS' stands for Year Start
plt.xlim(pd.Timestamp(start_date), pd.Timestamp(end_date))  # Set x-axis limits

plt.show()

