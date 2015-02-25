import centerlineMesh as FE
import os


FEmesh = FE.centerlineMesh(os.getcwd()+os.sep+'gmsh')
FEmesh.vtkReadVTKfile('SurfaceMesh.vtk')
FEmesh.getBoundary()
FEmesh.initializeSystem()
FEmesh.assembleK()
FEmesh.computeSol('art1.', 'art0.')
FEmesh.writeVTKsurface('SurfaceMesh')
FEmesh.computeCenterLine()
FEmesh.createVTKCenterLine('moy')
FEmesh.writeVTKCenterLine('CentreLine')
