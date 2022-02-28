# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 12:56:55 2022

@author: u0142106
"""
#DISCLAIMER!: THIS IS NOT A FINAL VERSION, PLEASE FEEL FREE TO SEND ME FEEDBACK 

#TODO: explore other combinations of available images from other satelites/spacecrafts
#So far C3 (easily modifiable), STA and STB Cor2. What happens if only 2 are available? is it possible to get other data
#PSP,SoLO, bepicolombo etc?
#It would be nice also if it could calculate speed (both total speed and radial speed) and time of insertion for spheromak

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
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.colors as colors


#Set your base and actual image to fit
#The base difference still is not incorporated, but it should be an easy task.
#Trying to see what is the best way of implementing
dateLASCO_base = dt.datetime(2013,3,12,12,18,0)
dateLASCO_img = dt.datetime(2013,3,12,12,45,0)

date_STA_base = dt.datetime(2013,3,12,7,54,0)
date_STA_img = dt.datetime(2013,3,12,12,54,0)

date_STB_base = dt.datetime(2013,3,12,7,54,0)
date_STB_img = dt.datetime(2013,3,12,12,54,0)


###################################################################
#Helioviewer image download
hv = helioviewer.HelioviewerClient()
for sourceid, obs in hv.data_sources.items():
    print(f"{sourceid}: {obs}") 

LASCO_0 = hv.download_jp2(dateLASCO_base, observatory= 'SOHO', instrument = 'LASCO', detector = 'C3')
LASCO_1 = hv.download_jp2(dateLASCO_img, observatory= 'SOHO', instrument = 'LASCO', detector = 'C3')

STA_0 = hv.download_jp2(date_STA_base, observatory='STEREO_A', instrument='SECCHI', detector='COR2')
STA_1 = hv.download_jp2(date_STA_img, observatory='STEREO_A', instrument='SECCHI', detector='COR2')

STB_0 = hv.download_jp2(date_STB_base,observatory='STEREO_B',instrument='SECCHI',detector='COR2')
STB_1 = hv.download_jp2(date_STB_img,observatory='STEREO_B',instrument='SECCHI',detector='COR2')
###################################################################

###################################################################
#processing images
###################################################################
cor2_A = sunpy.map.Map(STA_1)
cor2_B = sunpy.map.Map(STB_1)
C3 = sunpy.map.Map(LASCO_1)

sun_to_stereo_A = cor2_A.observer_coordinate.transform_to('hcrs')
sun_to_stereo_B = cor2_B.observer_coordinate.transform_to('hcrs')
sun_to_C3 = C3.observer_coordinate.transform_to('hcrs')

stereo_to_sun_A = SkyCoord(-sun_to_stereo_A.spherical, obstime=sun_to_stereo_A.obstime, frame='hcrs')
stereo_to_sun_B = SkyCoord(-sun_to_stereo_B.spherical, obstime=sun_to_stereo_B.obstime, frame='hcrs')
stereo_to_C3 = SkyCoord(-sun_to_C3.spherical, obstime=sun_to_C3.obstime, frame='hcrs')

###################################################################
#Maps difference
###################################################################
map_LASCO = Map([LASCO_0],[LASCO_1], sequence= True)
diff_LASCO = map_LASCO[1].data - map_LASCO[0].data
meta_LASCO = map_LASCO[1].meta
diff_map_LASCO = Map(diff_LASCO, meta_LASCO)

map_STA = Map([STA_0],[STA_1], sequence= True)
diff_STA = map_STA[1].data - map_STA[0].data
meta_STA = map_STA[1].meta
diff_map_STA = Map(diff_STA, meta_STA)

map_STB = Map([STB_0],[STB_1], sequence= True)
diff_STB = map_STB[1].data - map_STB[0].data
meta_STB = map_STB[1].meta
diff_map_STB = Map(diff_STB, meta_STB)
###################################################################


###################################################################
#PLANETS IN THE DIFFERENT FOVs
#TODO, add more planets and see if they should be plotted or not
mars_A = get_body_heliographic_stonyhurst('mars', cor2_A.date, observer=cor2_A.observer_coordinate)
mars_B = get_body_heliographic_stonyhurst('mars', cor2_B.date, observer=cor2_B.observer_coordinate)
mars_C3 = get_body_heliographic_stonyhurst('mars', C3.date, observer=C3.observer_coordinate)
###################################################################


def reset(event):
    Lon.reset()
    Lat.reset()
    Tilt.reset()
    Height.reset()
    Ratio.reset()
    Half_Angle.reset()
    straight_vertices.reset()
    front_vertices.reset()
    circle_vertices.reset()
    alpha.reset()
def update(val):
    #Updating parameters
    longitude_aux = (Lon.val/360)*2*np.pi 
    latitude_aux = (Lat.val/360)*2*np.pi
    Tilt_aux = (Tilt.val/360)*2*np.pi
    Height_aux = Height.val
    Ratio_aux = Ratio.val
    Half_Angle_aux = (Half_Angle.val/360)*2*np.pi
    straight_vertices_aux = straight_vertices.val
    front_vertices_aux = front_vertices.val
    circle_vertices_aux = circle_vertices.val
    alpha_aux = alpha.val
    marker_aux=marker.val
    textstr = '\n'.join((
    r'$LATITUDE = %.2f$' % (Lat.val, ),
    r'$LONGITUDE = %.2f$' % (Lon.val, ),
    r'$TILT = %.2f$' % (Tilt.val, ),
    r'$HEIGHT=%.2f$' % (Height.val, ),
   r'$RATIO (\delta) = %.2f$' % (Ratio.val, ),
   r'$\alpha = %.2f$' % (Half_Angle.val, ),
   r'$\delta = %.2f$' % (np.arcsin(Ratio.val), ),
   r'HALF WIDTH = %.2f' % (Half_Angle.val+np.arcsin(Ratio.val)*180/np.pi, ),))
    plt.text(0.55, -5.5,textstr,
         size=12,
         ha="right", bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),
                   )
         )
    print(longitude_aux,latitude_aux,Tilt_aux,Height_aux,Ratio_aux,Half_Angle_aux)
    ###################################################################
    #FIGURE GENERATOR
    #LASCO
    fig = plt.figure('C3',figsize= (5,5))
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
    mesh = gcs_mesh_sunpy(dateLASCO_img,
                          Half_Angle_aux,
                          Height_aux, 
                          straight_vertices_aux, 
                          front_vertices_aux,
                          circle_vertices_aux, 
                          Ratio_aux, 
                          latitude_aux, 
                          longitude_aux,
                          Tilt_aux)
    ax_2.plot_coord(mesh, '. ', ms=marker_aux,color = 'yellow',alpha = alpha_aux)
    fig.canvas.draw_idle()
    
    #STA
    fig = plt.figure('STA',figsize= (5,5))
    fig.clf()
    ax = plt.subplot(projection=cor2_A)
    # Let's tweak the axis to show in degrees instead of arcsec
    lon, lat = ax.coords
    lon.set_major_formatter('d.dd')
    lat.set_major_formatter('d.dd')
    cor2_A.plot(axes=ax, vmin=0, vmax=600)
    cor2_A.draw_limb()
    ax.plot_coord(mesh, '. ', ms=marker_aux,color = 'green',alpha = alpha_aux)
    # Plot the position of Mars
    #ax.plot_coord(mars_A, 's', color='white', fillstyle='none', markersize=12, label='Mars')
    # Plot all of the stars
    #ax.plot_coord(tbl_crds, 'o', color='white', fillstyle='none')
    fig.canvas.draw_idle()
       
    #STB
    fig = plt.figure('STB',figsize= (5,5))
    fig.clf()
    ax_1 = plt.subplot(projection=cor2_B)
    # Let's tweak the axis to show in degrees instead of arcsec
    lon, lat = ax_1.coords
    lon.set_major_formatter('d.dd')
    lat.set_major_formatter('d.dd')
    cor2_B.plot(axes=ax_1, vmin=0, vmax=600)
    cor2_B.draw_limb()
    ax_1.plot_coord(mesh, ' .', ms=marker_aux,color = 'green', alpha= alpha_aux)
    # Plot the position of Mars
    #ax_1.plot_coord(mars_B, 's', color='white', fillstyle='none', markersize=12, label='Mars')
    # Plot all of the stars
    #ax.plot_coord(tbl_crds, 'o', color='white', fillstyle='none')
    fig.canvas.draw_idle()
    plt.show()
###################################################################
#Figures without base difference
###################################################################

fig = plt.figure('Parameters',figsize=(5,10),tight_layout=True)
fig.clf()
#plt.subplots_adjust(left=0.25, bottom=0.25)
axcolor = 'lightgoldenrodyellow'
axlon = plt.axes([0.25, 0.95, 0.65, 0.03], facecolor=axcolor)
axlat = plt.axes([0.25, 0.90, 0.65, 0.03], facecolor=axcolor)
axtilt = plt.axes([0.25, 0.85, 0.65, 0.03], facecolor=axcolor)
axheight = plt.axes([0.25, 0.80, 0.65, 0.03], facecolor=axcolor)
axratio = plt.axes([0.25, 0.75, 0.65, 0.03], facecolor=axcolor)
axhalf = plt.axes([0.25, 0.70, 0.65, 0.03], facecolor=axcolor)
axstraight = plt.axes([0.25, 0.65, 0.65, 0.03], facecolor=axcolor)
axfront = plt.axes([0.25, 0.60, 0.65, 0.03], facecolor=axcolor)
axcircle = plt.axes([0.25, 0.55, 0.65, 0.03], facecolor=axcolor)
axalpha = plt.axes([0.25, 0.50, 0.65, 0.03], facecolor=axcolor)
axmarker = plt.axes([0.25, 0.45, 0.65, 0.03], facecolor=axcolor)


resetax = plt.axes([0.25, 0.40, 0.65, 0.03])
button = Button(resetax, 'Reset', color='red', hovercolor='0.975')
button.on_clicked(reset)

Lon = Slider(axlon, 'Longitude', 0, 359.9, valinit=0, valstep=0.1)
Lat = Slider(axlat, 'Latitude', -90, 90, valinit=0 , valstep= 0.1)
Tilt = Slider(axtilt, 'Tilt', -90, 90, valinit=0, valstep=0.1)
Height = Slider(axheight, 'Height', 0, 21.5, valinit=0.2 , valstep=0.1)
Ratio = Slider(axratio, 'Ratio', 0, 1, valinit=0.3)
Half_Angle = Slider(axhalf, 'Half Angle', 0, 90, valinit=40 , valstep=0.1)
straight_vertices = Slider(axstraight, 'Straight Vertices', 0, 100, valinit=5 , valstep=1)
front_vertices = Slider(axfront, 'Front Vertices', 0, 100, valinit=20 , valstep=1)
circle_vertices = Slider(axcircle, 'Circle Vertices', 0, 100, valinit=20 , valstep=1)
alpha = Slider(axalpha, 'Transparency', 0, 1, valinit=0.9 , valstep=0.01)
marker = Slider(axmarker, 'Marker Size', 1, 10, valinit=1 , valstep=1)
plt.tight_layout()


Lon.on_changed(update)
Lat.on_changed(update)
Tilt.on_changed(update)
Height.on_changed(update)
Ratio.on_changed(update)
Half_Angle.on_changed(update)
straight_vertices.on_changed(update)
front_vertices.on_changed(update)
circle_vertices.on_changed(update)
alpha.on_changed(update)
marker.on_changed(update)


plt.show()

