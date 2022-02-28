# GCS_tool GIT
Any comments, suggestions, ideas, modifications etc that you might have send them to me or edit the IDEAS.txt file (and notify me of a change as well please)
Hope this is useful!

## geometry.py:
This file is the one with the geometry for the GCS fit. 

## GCS_withSlider.py
This file is the one where you do the fitting.
Instrucitons: 
1) modify the datetime for the images you want to take and fit
2) run the program. A new window with the sliders should appear.
3) Modify something in the sliders and the images you selected should be prompted in the screen.

      -- at the begining the central part of the Sun should have a green/yellow mesh but very small. Modify the height and see how that works

4) Do the fit. The final parameters should appear in a red box below the reset button. 
5) Modify the transpaerncy and amount of points until you are happy.
6) Save the images at your deired location.

# IMPORTANT
The final result box is the one that has the parameters that should be input for CME models. Dont use the Half-Angle as the Half-Width because they are not the same!

## Image_viewer.py
This is a work in progress that would be useful for checking which images are available with the helioviewer tool for python before setting a date for doing the GCS fitting
