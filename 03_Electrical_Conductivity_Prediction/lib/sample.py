# Written by Marcus Schwarting in collaboration with the National Renewable
# Energy Laboratories. For further information please see the website 
# https://htem.nlr.gov/ as well as the recent publication on this work,
# [insert reference to upcoming paper here]. The goal of this project is to
# further the objectives laid out in the strategic objectives for the
# Materials Genome Initiative.
# 
# For further information on this code project, please email:
# Marcus Schwarting ==> marcus.schwarting@nrel.gov

#API POSITION(S) SPECIFIC!
import pandas as pd
from urllib.error import HTTPError
from library import get_htem_json

class Sample:
    
    def __init__(self,identity):
        self.identity = identity
        
    @staticmethod
    def search_by_ids(ids_list):
        obj_list = []
        for i in ids_list:
            obj_list.append(Sample(i))
        return obj_list

    def properties(self):
        try:
            data = get_htem_json('sample/' + str(self.identity))
        except HTTPError as error:
            if error.code in (400, 404):
                return pd.DataFrame()
            raise
        return pd.DataFrame([data])
            
    def spectra(self,which): 
        data = get_htem_json('sample/' + str(self.identity))
        #There is the potential to replace this with mvl_optical or mvl_xrd, 
        #but these seem to be broken at the moment...
        df = pd.DataFrame()
        
        leveled_position = 1
        if which == 'xrf':
            try:
                #b = len(df)
                #for i in range(len(data['xrf_compounds'])):
                #    df.at[b,'xrf_concentrations_'+(data['xrf_compounds'][i])] = data['xrf_concentration'][i]     
                df['xrf_compounds'] = data['xrf_compounds']
                df['xrf_concentrations'] = data['xrf_concentration']
            except:
                #display('No xrf data')
                pass
        elif which == 'xrd': 
            try:
                df['xrd_angle'] = data['xrd_angle']
                df['xrd_background'] = data['xrd_background']
                df['xrd_intensity'] = data['xrd_intensity']
            except:
                #display('No xrd data')
                pass
        elif which == 'opt':
            info_df = pd.DataFrame()
            extra_df = pd.DataFrame()
            try:
                info_df['peak_count'] = list([data['peak_count']])
                info_df['opt_average_vis_trans'] = list([data['opt_average_vis_trans']])
                info_df['opt_direct_bandgap'] = list([data['opt_direct_bandgap']])
                
                extra_df['xyz_mm'] = data['xyz_mm']
            except:
                pass
            uvir_df = pd.DataFrame()
            try:
                uvir_df['opt_uvir_wavelength'] = data['opt_uvir_wavelength']
                uvir_df['opt_uvir_response'] = data['opt_uvir_response']
            except KeyError: #No uvir available
                pass
            uvit_df = pd.DataFrame()
            try:
                uvit_df['opt_uvit_wave'] = data['opt_uvit_wavelength']
                uvit_df['opt_uvit_response'] = data['opt_uvit_response']
            except KeyError: #No uvit available
                pass
            nirr_df = pd.DataFrame()
            try:
                nirr_df['opt_nirr_wavelength'] = data['opt_nirr_wavelength']
                nirr_df['opt_nirr_response'] = data['opt_nirr_response']
            except KeyError: #No nirr available
                pass
            nirt_df = pd.DataFrame()
            try:
                nirt_df['opt_nirt_wavelength'] = data['opt_nirt_wavelength']
                nirt_df['opt_nirt_response'] = data['opt_nirt_response']
            except KeyError: #No nirt available
                pass
            opt_df = pd.DataFrame()
            try:
                opt_df['opt_energy'] = data['opt_energy']
                opt_df['opt_wavelength'] = data['opt_wavelength']
                opt_df['opt_normalized_transmittance'] = data['opt_normalized_transmittance']
                opt_df['opt_absorption_coefficient'] = data['opt_absorption_coefficient']
            except: #No opt available
                #display('No optical data')
                pass
            pos_df = pd.concat([info_df,extra_df,uvir_df,uvit_df,nirr_df,nirt_df,opt_df],axis=1) 
            df = pd.concat([df,pos_df],axis=1)
        elif which == 'ele':
            ele_df = pd.DataFrame()
            info_df = pd.DataFrame()
            try:
                ele_df['fpm_voltage_volts'] = data['fpm_voltage_volts'] #5
                ele_df['fpm_current_amps'] = data['fpm_current_amps'] #5
                
                info_df['fpm_sheet_resistance'] = list([data['fpm_sheet_resistance']]) #1 
                info_df['fpm_standard_deviation'] = list([data['fpm_standard_deviation']]) #1
                info_df['fpm_resistivity'] = list([data['fpm_resistivity']]) #1
                info_df['fpm_conductivity'] = list([data['fpm_conductivity']]) #1
                info_df['absolute_temp_c'] = list([data['absolute_temp_c']]) 
                #if data['STRING'] == int else 
            except:
                #display('No electrical data')
                pass
            pos_df = pd.concat([ele_df,info_df],axis=1)
            df = pd.concat([df,pos_df],axis=1)
        else:
            df = pd.DataFrame()
        return df
