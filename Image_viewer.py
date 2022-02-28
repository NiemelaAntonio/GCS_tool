# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 16:03:45 2022

@author: u0142106
"""
import numpy as np
import matplotlib.pyplot as plt
import datetime
from datetime import datetime
import sys
#sys.path.append('C:/Users/u0142106/Desktop/New folder/GCS_python/')
from geometry import gcs_mesh_sunpy
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import Distance, SkyCoord
from astropy.time import Time
import sunpy.map
from sunpy.map import Map
from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.net import helioviewer
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import matplotlib.colors as colors

hv = helioviewer.HelioviewerClient()

    
Year_=['2013']
Month_=['02']
Day_=['01']
Hour_=['00']
Minute_=['00']
Second_=['00']   
datetime_input=[]
datetime_input_aux=[]
spacecraft=[]


class DATE:

    def Year(text):
        Year_.clear()
        datetime_input.clear()
        Year_.append(text)
        datetime_input.append(str(Year_[0])+str('-')+str(Month_[0])+str('-')+str(Day_[0])+str('T')
                        +str(Hour_[0])+str(':')+str(Minute_[0])+str(':')+str(Second_[0]))
        print(datetime_input)
    def Month(text):
        Month_.clear()
        datetime_input.clear()
        Month_.append(text)
        datetime_input.append(str(Year_[0])+str('-')+str(Month_[0])+str('-')+str(Day_[0])+str('T')
                        +str(Hour_[0])+str(':')+str(Minute_[0])+str(':')+str(Second_[0]))
        print(datetime_input)
    def Day(text):
        Day_.clear()
        datetime_input.clear()
        Day_.append(text)
        datetime_input.append(str(Year_[0])+str('-')+str(Month_[0])+str('-')+str(Day_[0])+str('T')
                        +str(Hour_[0])+str(':')+str(Minute_[0])+str(':')+str(Second_[0]))
        print(datetime_input)
    def Hour(text):
        Hour_.clear()
        datetime_input.clear()
        Hour_.append(text)
        datetime_input.append(str(Year_[0])+str('-')+str(Month_[0])+str('-')+str(Day_[0])+str('T')
                        +str(Hour_[0])+str(':')+str(Minute_[0])+str(':')+str(Second_[0]))
        print(datetime_input)
    def Minute(text):
        Minute_.clear()
        datetime_input.clear()
        Minute_.append(text)
        datetime_input.append(str(Year_[0])+str('-')+str(Month_[0])+str('-')+str(Day_[0])+str('T')
                        +str(Hour_[0])+str(':')+str(Minute_[0])+str(':')+str(Second_[0]))
        print(datetime_input)
    def Second(text):
        Second_.clear()
        Second_.append(text)
        datetime_input.clear()
        datetime_input.append(str(Year_[0])+str('-')+str(Month_[0])+str('-')+str(Day_[0])+str('T')
                        +str(Hour_[0])+str(':')+str(Minute_[0])+str(':')+str(Second_[0]))
        print(datetime_input)
    datetime_input_aux=datetime_input
    print(datetime_input)
        #print(type(datetime_input[0]))
        #print(datetime_input[0])


       
class IMAGES:
    spacecraft=[]
    def C3(event):
        spacecraft.clear()
        spacecraft.append('C3')
        print(spacecraft)
        
    def C2(event):
        spacecraft.clear()
        spacecraft.append('C2')
        print(spacecraft)
        
    def STA(event):
        spacecraft.clear()
        spacecraft.append('STA')
        print(spacecraft)
                
    def STB(event):
        spacecraft.clear()
        spacecraft.append('STB')
        print(spacecraft)
    def plotting(event):
        print(datetime_input[0])
        datetime_image=datetime.strptime(datetime_input[0],'%Y-%m-%dT%H:%M:%S')
        if spacecraft[0]== 'C3':
            print('getting image')
            LASCO_C3 = hv.download_jp2(datetime_image, observatory= 'SOHO', instrument = 'LASCO', detector = 'C3')
            if len(LASCO_C3)>0:
                print('image saved in '+str(LASCO_C3))
            else:
                print('Nothing  found')
            C3 = sunpy.map.Map(LASCO_C3)
            sun_to_C3 = C3.observer_coordinate.transform_to('hcrs')
            stereo_to_C3 = SkyCoord(-sun_to_C3.spherical, obstime=sun_to_C3.obstime, frame='hcrs')
            print('Image coming out of the oven')
            fig = plt.figure('C3',figsize= (7.5,7.5))
            fig.clf()
            ax_2 = plt.subplot(projection=C3)
            # Let's tweak the axis to show in degrees instead of arcsec
            lon, lat = ax_2.coords
            lon.set_major_formatter('d.dd')
            lat.set_major_formatter('d.dd')
            C3.plot_settings['cmap'] = plt.get_cmap('Blues_r')
            C3.plot_settings['norm'] = colors.LogNorm(100, C3.max())
            C3.plot(axes=ax_2)
            C3.draw_limb()
            fig.canvas.draw_idle()
            plt.show()
            
            
        if spacecraft[0]== 'C2':
            LASCO_C2 = hv.download_jp2(datetime_image, observatory= 'SOHO', instrument = 'LASCO', detector = 'C2')
            print('getting image')
            if len(LASCO_C2)>0:
                print('image saved in '+str(LASCO_C2))
            else:
                print('Nothing  found')
            C2 = sunpy.map.Map(LASCO_C2)
            sun_to_C2 = C2.observer_coordinate.transform_to('hcrs')
            stereo_to_C2 = SkyCoord(-sun_to_C2.spherical, obstime=sun_to_C2.obstime, frame='hcrs')
            fig = plt.figure('C2',figsize= (7.5,7.5))
            fig.clf()
            ax_2 = plt.subplot(projection=C2)
            # Let's tweak the axis to show in degrees instead of arcsec
            lon, lat = ax_2.coords
            lon.set_major_formatter('d.dd')
            lat.set_major_formatter('d.dd')
            C2.plot_settings['cmap'] = plt.get_cmap('Oranges_r')
            C2.plot_settings['norm'] = colors.LogNorm(100, C2.max())
            C2.plot(axes=ax_2)
            C2.draw_limb()
            fig.canvas.draw_idle()
            plt.show()
            
            
        if spacecraft[0]== 'STA':
            print('getting image')
            STA_ = hv.download_jp2(datetime_image, observatory='STEREO_A', instrument='SECCHI', detector='COR2')
            print('getting image')
            if len(STA_)>0:
                print('image saved in '+str(STA_))
            else:
                print('Nothing  found')
            cor2_A = sunpy.map.Map(STA_)
            sun_to_stereo_A = cor2_A.observer_coordinate.transform_to('hcrs')
            stereo_to_sun_A = SkyCoord(-sun_to_stereo_A.spherical, obstime=sun_to_stereo_A.obstime, frame='hcrs')
            fig = plt.figure('STA',figsize= (7.5,7.5))
            fig.clf()
            ax = plt.subplot(projection=cor2_A)
            # Let's tweak the axis to show in degrees instead of arcsec
            lon, lat = ax.coords
            lon.set_major_formatter('d.dd')
            lat.set_major_formatter('d.dd')
            cor2_A.plot(axes=ax, vmin=0, vmax=600)
            cor2_A.draw_limb()
            fig.canvas.draw_idle()
            plt.show()
    
            
        if spacecraft[0]== 'STB':
            STB_ = hv.download_jp2(datetime_image, observatory='STEREO_B', instrument='SECCHI', detector='COR2')
            print('getting image')
            if len(STB_)>0:
                print('image saved in '+str(STB_))
            else:
                print('Nothing  found')
            cor2_B = sunpy.map.Map(STB_)
            sun_to_stereo_B = cor2_B.observer_coordinate.transform_to('hcrs')
            stereo_to_sun_B = SkyCoord(-sun_to_stereo_B.spherical, obstime=sun_to_stereo_B.obstime, frame='hcrs')
            fig = plt.figure('STB',figsize= (7.5,7.5))
            fig.clf()
            ax_1 = plt.subplot(projection=cor2_B)
            # Let's tweak the axis to show in degrees instead of arcsec
            lon, lat = ax_1.coords
            lon.set_major_formatter('d.dd')
            lat.set_major_formatter('d.dd')
            cor2_B.plot(axes=ax_1, vmin=0, vmax=600)
            cor2_B.draw_limb()
            fig.canvas.draw_idle()
            plt.show()

def reset(event):
    spacecraft.clear()
    datetime_input.clear()




ax = plt.figure('Choosing Images',figsize=(10,5))
#plt.subplots_adjust(bottom=0.2)
initial_text = '%Y-%m-%dT%H:%M:%S'
    
axC3 = plt.axes([0.1, 0.85, 0.1, 0.075])
axC2 = plt.axes([0.21, 0.85, 0.1, 0.075])
axSTA = plt.axes([0.32, 0.85, 0.1, 0.075])
axSTB = plt.axes([0.43, 0.85, 0.1, 0.075])
bC3 = Button(axC3, 'C3',hovercolor='0.95',color='grey')
bC3.on_clicked(IMAGES.C3)
bC2 = Button(axC2, 'C2',hovercolor='0.95',color='grey')
bC2.on_clicked(IMAGES.C2)
bSTA = Button(axSTA, 'STA',hovercolor='0.95',color='grey')
bSTA.on_clicked(IMAGES.STA)
bSTB = Button(axSTB, 'STB',hovercolor='0.95',color='grey')
bSTB.on_clicked(IMAGES.STB)

axbox_year = plt.axes([0.33, 0.75, 0.1, 0.05])
text_box_1 = TextBox(axbox_year, 'Year (4 digits)', initial='2013')
text_box_1.on_submit(DATE.Year)
axbox_month = plt.axes([0.33, 0.65, 0.1, 0.05])
text_box_2 = TextBox(axbox_month, 'Month (2 digits)', initial='05')
text_box_2.on_submit(DATE.Month)
axbox_day = plt.axes([0.33, 0.55, 0.1, 0.05])
text_box_3 = TextBox(axbox_day, 'Day (2 digits)', initial='01')
text_box_3.on_submit(DATE.Day)
axbox_hour = plt.axes([0.33, 0.45, 0.1, 0.05])
text_box_4 = TextBox(axbox_hour, 'Hour (2 digits)', initial='00')
text_box_4.on_submit(DATE.Hour)
axbox_minute = plt.axes([0.33, 0.35, 0.1, 0.05])
text_box_5 = TextBox(axbox_minute, 'Minute (2 digits)', initial='00')
text_box_5.on_submit(DATE.Minute)
axbox_second = plt.axes([0.33, 0.25, 0.1, 0.05])
text_box_6 = TextBox(axbox_second, 'Seconds (2 digits)', initial='00')
text_box_6.on_submit(DATE.Second)



resetax = plt.axes([0.15, 0.1, 0.1, 0.05])
button = Button(resetax, 'Reset', color='red', hovercolor='0.975')
button.on_clicked(reset)

Plotting = plt.axes([0.30, 0.1, 0.1, 0.05])
button_1 = Button(Plotting, 'Plot!', color='red', hovercolor='0.975')
button_1.on_clicked(IMAGES.plotting)

