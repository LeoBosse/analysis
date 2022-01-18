import math
import time
import vtk
import numpy as np
import matplotlib.pyplot as plt



def getPointH(A,r,altitude,elevation,azimuth):
   ele=math.radians(elevation)
   azi=math.radians(azimuth)
   h=altitude
   lonA=math.radians(A[0])
   latA=math.radians(A[1])
   d=-r*math.sin(ele)+math.sqrt((r**2)*(math.sin(ele)**2)+h**2+2*r*h)
   beta=math.atan2(d*math.cos(ele),r+d*math.sin(ele))
   # latH and lonH calculations taken from
   #  https://www.movable-type.co.uk/scripts/latlong.html
   latH = math.asin(math.sin(latA)*math.cos(beta) +
          math.cos(latA)*math.sin(beta)*math.cos(azi))
   lonH = lonA + math.atan2(math.sin(azi)*math.sin(beta)*math.cos(latA),
                            math.cos(beta)-math.sin(latA)*math.sin(latH))
   H = [math.degrees(lonH),math.degrees(latH)]
   return H

def toCartesian(r,lon,lat):
   theta=math.radians(90.0-lat)
   phi=math.radians(lon)
   x=r*math.sin(theta)*math.cos(phi)
   y=r*math.sin(theta)*math.sin(phi)
   z=r*math.cos(theta)
   return x,y,z

# Earth radius
r=6371.0
# Ionospheric observation point altitude
altitude = 220.0

# Create a sphere Earth
sphereSource = vtk.vtkSphereSource()
sphereSource.SetCenter(0.0, 0.0, 0.0)
sphereSource.SetRadius(1.0)
sphereSource.SetThetaResolution(90)
sphereSource.SetPhiResolution(45)
sphereSource.Update()
mapper = vtk.vtk.vtkPolyDataMapper()
mapper.SetInputConnection(sphereSource.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
#actor.GetProperty().SetRepresentationToWireframe()
actor.GetProperty().SetColor(0.0, 1.0, 0.0)


# Create a sphere ionosphere
sphere2Source = vtk.vtkSphereSource()
sphere2Source.SetCenter(0.0, 0.0, 0.0)
sphere2Source.SetRadius(1.0+(altitude/r))
sphere2Source.SetThetaResolution(90)
sphere2Source.SetPhiResolution(45)
sphere2Source.Update()
mapper2 = vtk.vtk.vtkPolyDataMapper()
mapper2.SetInputConnection(sphere2Source.GetOutputPort())
actor2 = vtk.vtkActor()
actor2.SetMapper(mapper2)
actor2.GetProperty().SetRepresentationToWireframe()
actor2.GetProperty().SetColor(0.0, 0.0, 1.0)

# Place camera
camera = vtk.vtkCamera()


# Create a renderer, render window, and interactor
renderer = vtk.vtkRenderer()
renderer.SetActiveCamera(camera)

renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindow.SetSize(1280, 800)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderWindowInteractor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()


# Point A ==> Location of the SPP camera
lonA = 11.63
latA = 78.62
A = [lonA,latA]
print( "A degrees: ",A)
Ac=toCartesian(r,A[0],A[1])
Ac_norm=[Ac[0]/r,Ac[1]/r,Ac[2]/r,]


n_steps = 90
H_lons = np.zeros((4, n_steps+1), dtype=float)
H_lats = np.zeros((4, n_steps+1), dtype=float)
AoLP = np.zeros((4, n_steps+1), dtype=float)
elevs = np.arange(n_steps+1)

azis = np.arange(4)
for a in azis:

        # Camera direction
        azimuth = a*90.0

        for e in elevs:
           el=e*1.
           print( "Elevation: ",el   )
           print( "Azimuth: ", azimuth  )
           # Point H ==> Observation point
           H = getPointH(A,r,altitude,el,azimuth)
           H_lons[a,e] = H[0]
           H_lats[a,e] = H[1]
           print( H[0], H[1])

           Hc=toCartesian(r+altitude,H[0],H[1])
           Hc_norm=[Hc[0]/r,Hc[1]/r,Hc[2]/r,]


           # Create the line of sight of the camera
           line=vtk.vtkLineSource()
           line.SetPoint1(Ac_norm)
           line.SetPoint2(Hc_norm)
           linemapper = vtk.vtk.vtkPolyDataMapper()
           linemapper.SetInputConnection(line.GetOutputPort())
           lineactor = vtk.vtkActor()
           lineactor.SetMapper(linemapper)
           lineactor.GetProperty().SetColor(1.0, 0.0+el/80., 0.0)

           # Create a line representing magnetic field
           magpoint1=Hc_norm
           magpoint2=[Hc_norm[0]+0.04,Hc_norm[1]+0.04,Hc_norm[2]+0.04]
           arrow=vtk.vtkLineSource()
           arrow.SetPoint1(magpoint1)
           arrow.SetPoint2(magpoint2)
           arrowmapper = vtk.vtk.vtkPolyDataMapper()
           arrowmapper.SetInputConnection(arrow.GetOutputPort())
           arrowactor = vtk.vtkActor()
           arrowactor.SetMapper(arrowmapper)
           arrowactor.GetProperty().SetColor(0.0, 0.0, 1.0)

           # Change camera
           camera.SetPosition(Ac_norm[0],Ac_norm[1],Ac_norm[2])
           #camera.Azimuth(azimuth)
           #camera.Elevation(el)
           camera.SetFocalPoint(Hc_norm[0],Hc_norm[1],Hc_norm[2])
           camera.SetViewUp(Ac_norm[0],Ac_norm[1],Ac_norm[2])
           camera.OrthogonalizeViewUp

           # Add the actors to the scene
           renderer.AddActor(lineactor)
           renderer.AddActor(arrowactor)


           # Calculate view coordinates
           coor1=vtk.vtkCoordinate()
           coor2=vtk.vtkCoordinate()
           coor1.SetValue(magpoint1[0],magpoint1[1],magpoint1[2])
           coor2.SetValue(magpoint2[0],magpoint2[1],magpoint2[2])
           mag1view=coor1.GetComputedViewportValue(renderer)
           mag2view=coor2.GetComputedViewportValue(renderer)
           #AoLP[a,e]= math.degrees(math.atan2(mag2view[1]-mag1view[1],mag2view[0]-mag1view[0]))
           AoLP[a,e]= math.degrees(math.atan2(mag2view[0],mag2view[1]))
           print( AoLP[a,e])
           renderWindow.Render()
           # Add the sphere to the scene
           renderer.AddActor(actor)
           renderer.AddActor(actor2)
           renderWindowInteractor.Start()



# Plot the results
plt.subplot(3,1,1)
N_azis_lon = H_lons[0,:]
E_azis_lon = H_lons[1,:]
S_azis_lon = H_lons[2,:]
W_azis_lon = H_lons[3,:]
plt.plot(elevs, N_azis_lon, 'r--', elevs, E_azis_lon, 'g.', elevs, S_azis_lon, 'b^', elevs, W_azis_lon, 'ys' )


plt.subplot(3,1,2)
N_azis_lat = H_lats[0,:]
E_azis_lat = H_lats[1,:]
S_azis_lat = H_lats[2,:]
W_azis_lat = H_lats[3,:]
plt.plot(elevs, N_azis_lat, 'r--', elevs, E_azis_lat, 'g.', elevs, S_azis_lat, 'b^', elevs, W_azis_lat, 'ys' )

plt.subplot(3,1,3)
plt.plot(elevs, AoLP[0,:], 'r--', elevs, AoLP[1,:], 'g.', elevs, AoLP[2,:], 'b^', elevs, AoLP[3,:], 'ys' )

plt.show()

# Add the sphere to the scene
renderer.AddActor(actor)
renderer.AddActor(actor2)

# Set final camera position
#camera.SetPosition(0,-4,4)
#camera.SetFocalPoint(Ac_norm[0],Ac_norm[1],Ac_norm[2],)
#camera.SetViewUp(0,0,1)

# Render and interact
renderWindow.Render()



renderWindowInteractor.Start()
