# -*- coding: utf-8 -*-
#!===========================================================================
#!
#! Copyright 2013-2015 CNRS.
#!
#! This file is part of a software (dicom2stl) which is a computer program 
#! which purpose is to modelize mesh from DICOM files.
#!
#! This software is governed by the CeCILL license under French law and
#! abiding by the rules of distribution of free software.  You can  use, 
#! modify and/ or redistribute the software under the terms of the CeCILL
#! license as circulated by CEA, CNRS and INRIA at the following URL
#! "http://www.cecill.info". 
#!
#! As a counterpart to the access to the source code and  rights to copy,
#! modify and redistribute granted by the license, users are provided only
#! with a limited warranty  and the software's author,  the holder of the
#! economic rights,  and the successive licensors  have only  limited
#! liability. 
#!
#! In this respect, the user's attention is drawn to the risks associated
#! with loading,  using,  modifying and/or developing or reproducing the
#! software by the user in light of its specific status of free software,
#! that may mean  that it is complicated to manipulate,  and  that  also
#! therefore means  that it is reserved for developers  and  experienced
#! professionals having in-depth computer knowledge. Users are therefore
#! encouraged to load and test the software's suitability as regards their
#! requirements in conditions enablinfrom vtk.util import numpy_support as VNg the security of their systems and/or 
#! data to be ensured and,  more generally, to use and operate it in the 
#! same conditions as regards security. 
#!
#! The fact that you are presently reading this means that you have had
#! knowledge of the CeCILL license and that you accept its terms.
#!
#! To report bugs, suggest enhancements, etc. to the Authors, contact
#! Dominique Ambard.
#!
#! dominique.ambard@um2.fr
#!
#!===========================================================================
import numpy as np
import vtk,os
import scipy.sparse as sp
from vtk.util import numpy_support as VN 
#----------------------------------------------------------------------------------------------  
class centerlineMesh():
   """ Classe de fonction pour realiser des calculs de thermique staitonnaire """
   def __init__(self, read_path, log):
      """ Initialisation de la classe """
      self.read_path = read_path
      self.log = log
      self.deriveFonction()
      self.matConduction()
      self.gaussPoint()
      self.resolution = 100
      self.filter = 5
      self.listBoundary = []
      
   def vtkReadVTKfile(self, filename = 'not defined'):
      """ Lecture d'un fichier vtk"""
      # Read the source file.
      print "----------------------------------------------------"
      print "Read geometric surface"
      print "----------------------------------------------------"
      reader = vtk.vtkUnstructuredGridReader()
      reader.SetFileName(os.path.join(self.read_path,filename))
      reader.Update() # Needed because of GetScalarRange
      
      surface_filter = vtk.vtkDataSetSurfaceFilter()
      surface_filter.SetInputConnection(reader.GetOutputPort()) 
      
      self.meshNormals = vtk.vtkPolyDataNormals()
      self.meshNormals.SetInputConnection(surface_filter.GetOutputPort())
      self.meshNormals.Update()
      
      #~ self.normales = VN.vtk_to_numpy(self.meshNormals.GetOutput().GetPointData().GetNormals())
      self.meshNode  = VN.vtk_to_numpy(surface_filter.GetOutput().GetPoints().GetData())
      self.meshTable = VN.vtk_to_numpy(surface_filter.GetOutput().GetPolys().GetData())
      self.meshTable = self.meshTable.reshape(self.meshTable.shape[0]/4, 4)[:,1:]
      
      self.nbDof = self.meshNode.shape[0]
      self.nbElem = self.meshTable.shape[0]
      
   def setMesh(self, nodes, elements):
      "Affecte un maillage"
      
      print "----------------------------------------------------"
      print "Apply surface"
      print "----------------------------------------------------"
      
      self.meshNode = nodes
      self.meshTable = elements.reshape(elements.shape[0]/4, 4)[:,1:]
      self.nbDof = self.meshNode.shape[0]
      self.nbElem = self.meshTable.shape[0]
      
      np.save('connec.npy',self.meshTable)
      
      print "----------------------------------------------------"
      print " nb nodes   : ",self.nbDof
      print " nb bulks   : ",self.nbElem
      print "----------------------------------------------------"
      
   def getBoundary(self):
      "Recherche des frontieres"
      print "----------------------------------------------------"
      print "Get boundary "
      print "----------------------------------------------------"
      # Recherche des frontieres dans le domaine
      boundaryEdge = []
      for i in range(self.meshTable.shape[0]):
         tri = self.meshTable[i,:]
         if len((np.sum((self.meshTable==tri[0])*1 + (self.meshTable==tri[1])*1,axis = 1) == 2).nonzero()[0])==1 :
            boundaryEdge.append([tri[0], tri[1]])
         if len((np.sum((self.meshTable==tri[1])*1 + (self.meshTable==tri[2])*1,axis = 1) == 2).nonzero()[0])==1 :
            boundaryEdge.append([tri[1], tri[2]])
         if len((np.sum((self.meshTable==tri[2])*1 + (self.meshTable==tri[0])*1,axis = 1) == 2).nonzero()[0])==1 :
            boundaryEdge.append([tri[2], tri[0]])
      
      boundaryEdge = np.array(boundaryEdge)
      
      # Decomposition des frontieres
      boundary = {}
      j = 0
      while len(boundaryEdge)>0 :
         boundary['art'+str(j)+'.'] = []
         
         boundary['art'+str(j)+'.'].append(boundaryEdge[0,0])
         boundary['art'+str(j)+'.'].append(boundaryEdge[0,1])
         
         node = boundaryEdge[0,1]
         boundaryEdge = np.delete(boundaryEdge,0,axis = 0)
         while node <> boundary['art'+str(j)+'.'][0]:
            lig = (boundaryEdge == node).nonzero()[0][0]
            col = (boundaryEdge == node).nonzero()[1][0]
            if col == 0 :
               boundary['art'+str(j)+'.'].append(boundaryEdge[lig,1])
               node = boundaryEdge[lig,1]
            if col == 1 :
               boundary['art'+str(j)+'.'].append(boundaryEdge[lig,0])
               node = boundaryEdge[lig,0]
            boundaryEdge = np.delete(boundaryEdge,lig,axis = 0)
         boundary['art'+str(j)+'.'].remove(boundary['art'+str(j)+'.'][0])
         j += 1
         
      self.boundarys = boundary
      
      print "----------------------------------------------------"
      print "Some boundary are detect"
      for name in self.boundarys.keys():
         print "Boundary - ",name
         #~ self.createVTKContour( self.meshNode[self.boundarys[name],:] )
         #~ self.writeVTKContour( name+'vtp')
         self.listBoundary.append(name)
      print "----------------------------------------------------"
      
   def formeFonction(self,x,y):
      "Fonction de forme d'un triangle 3 noeuds"
      
      self.N = np.zeros(3,float)
      self.N[0] = 1.0 - x - y
      self.N[1] = x
      self.N[2] = y
      
   def deriveFonction(self):
      "Derive des fonctions de formes"
      
      self.dN =np.zeros((3,2),float)
      self.dN[0,:] = np.array([-1.0,-1.0])
      self.dN[1,:] = np.array([ 1.0, 0.0])
      self.dN[2,:] = np.array([ 0.0, 1.0])
      
   def gaussPoint(self):
      "Applique la methode d'integration"
      
      self.PG = np.array([1.0/3.0, 1.0/3.0, 1.0/2.0])
   
   def matPassageLocal2Global(self,X):
      "Determination de la matrice de passage"
      
      Xcenter = np.mean(X,axis = 0)
      X = X - Xcenter
      x1 = X[1,:] - X[0,:]
      x1 = x1 / np.linalg.norm(x1)
      x2 = X[2,:] - X[0,:]
      x2 = x2 / np.linalg.norm(x2)
      n  = np.cross(x1,x2)
      n = n / np.linalg.norm(n)
      s  = x1
      t  = np.cross(n,s)
      
      P_LG = np.zeros((3,3),float)
      P_LG[:,0] = s
      P_LG[:,1] = t
      P_LG[:,2] = n

      return P_LG
   
   def matConduction(self):
      "Calcule la matrice de conduction"
      
      self.D = np.array([[1.0,0.0],[0.0,1.0]])
      
   def gradientIso(self,X): 
      "Calcule la matrice gradient"
      
      J = np.zeros((2,2),float)
      
      J[0,0] = np.dot(self.dN[:,0],X[:,0])
      J[0,1] = np.dot(self.dN[:,0],X[:,1])
      
      J[1,0] = np.dot(self.dN[:,1],X[:,0])
      J[1,1] = np.dot(self.dN[:,1],X[:,1])
      
      detJ = np.linalg.det(J)
      invJ = np.linalg.inv(J)
      
      B = np.dot(invJ,self.dN.T)
      
      return B,detJ
      
   def matElem(self,X):
      "Calcul de la matrice elementaire"
      
      B,detJ = self.gradientIso(X)
      Ke = self.PG[2]*detJ*(np.dot(B.T,np.dot(self.D,B)))
      
      return Ke
      
   def getNormFlux(self,X,dof):
      "Calcul les flux"
      
      B,detJ = self.gradientIso(X)
      Flux = np.linalg.norm(np.dot(B,dof))
      
      return Flux
   
   def initializeSystem(self):
      "Initialise le systeme matriciel"
      
      # Initialisation du systeme matriciel
      self.A = np.zeros((self.nbDof, self.nbDof), float)
      self.X = np.zeros(self.nbDof,float)
      self.Y = np.zeros(self.nbDof,float)

   def assembleK(self):
      "Creation de la matrice totale"
      print "----------------------------------------------------"
      print "Compute system"
      print "----------------------------------------------------"
      for i in range(self.nbElem):
         X = self.meshNode[self.meshTable[i,:],:]
         P = self.matPassageLocal2Global(X)
         Xloc = 0.0*X
         Xloc[0,:] = np.dot(P.T,X[0,:])
         Xloc[1,:] = np.dot(P.T,X[1,:])
         Xloc[2,:] = np.dot(P.T,X[2,:])
         Ke = self.matElem(Xloc[:,:2])
         
         for j,ind in enumerate(self.meshTable[i,:]):
            self.A[self.meshTable[i,:],ind] += Ke[j,:]
   
   def postFlux(self):
      "post traitement des flux nodaux"
      print "----------------------------------------------------"
      print "Post system"
      print "----------------------------------------------------"
      # boucle sur les elements
      Q = np.zeros(self.nbElem,float)
      for i in range(self.nbElem):
         X = self.meshNode[self.meshTable[i,:],:]
         P = self.matPassageLocal2Global(X)
         Xloc = 0.0*X
         Xloc[0,:] = np.dot(P.T,X[0,:])
         Xloc[1,:] = np.dot(P.T,X[1,:])
         Xloc[2,:] = np.dot(P.T,X[2,:])
         Q[i] = self.getNormFlux(Xloc[:,:2],self.X[self.meshTable[i,:]])
      
      # Boucle sur les noeuds
      for i in range(self.nbDof):
         self.Y[i] = np.mean(Q[self.getElemfromNode(i)])
   
   def getElemfromNode(self, node):
      "recupere les elements attenant a un noeud"
      return(self.meshTable == node).nonzero()[0]
   
   def computeSol(self, bd1, bd2):
      "Calcule la solution entre bd1 et bd2"
      
      n1 = self.boundarys[bd1]
      n2 = self.boundarys[bd2]
      
      # definition de dof connu
      connu = []
      connu.extend(n1)
      connu.extend(n2)
      
      # definition de dof inconnu
      inconnu = range(self.nbDof)
      for i in connu:
         try :
            inconnu.remove(i)
         except:
            pass
            
      # recuperation des matrices
      Kii = self.A[:,inconnu]
      Kii = Kii[inconnu,:]
      Kic = self.A[inconnu,:]
      Kic = Kic[:,connu]
      
      # determination des temperatures imposees
      Tc = np.concatenate((np.ones(len(n1),float),-np.ones(len(n2),float)),axis = 0)
      Q = - np.dot(Kic,Tc)
      
      # resolution du probleme
      Ti = np.linalg.solve(Kii, Q)
      
      # affectation des resultats
      for i,ind in enumerate(inconnu):
         self.X[ind] = Ti[i]
      for i,ind in enumerate(connu):
         self.X[ind] = Tc[i]
         
      # post traitement des flux
      self.postFlux()
      
   def writeVTKsurface(self, filename):
      "Create data in VTK format to compute isocontour"
      
      # Conteneur du nouveau maillage
      data = vtk.vtkUnstructuredGrid()
      pts = vtk.vtkPoints()
      for i in range(self.nbDof):
         pts.InsertPoint(i,self.meshNode[i,0], self.meshNode[i,1], self.meshNode[i,2])
      data.SetPoints(pts)
      for i in range(self.nbElem):
         Tri3 =vtk.vtkTriangle()
         for j,inode in enumerate(self.meshTable[i,:]):
            Tri3.GetPointIds().SetId(j, inode)
         data.InsertNextCell(Tri3.GetCellType(),Tri3.GetPointIds())
      
      value0 = vtk.vtkDoubleArray()
      value0.SetName("Value0")
      value0.SetNumberOfComponents(1)
      for i in range(self.nbDof):
         value0.InsertTuple1(i, self.X[i])
      data.GetPointData().AddArray(value0)
      
      value1 = vtk.vtkDoubleArray()
      value1.SetName("Value1")
      value1.SetNumberOfComponents(1)
      for i in range(self.nbDof):
         value1.InsertTuple1(i, self.Y[i])
      data.GetPointData().AddArray(value1)
      
      # Ecriture du nouveau maillage
      writer = vtk.vtkXMLUnstructuredGridWriter()
      writer.SetInput(data)
      writer.SetFileName(filename+'.vtu')
      writer.Update()
      
   def createVTKpolydata(self, val):
      "Create polydata VTK format to compute isocontour"
      
      vtkcoords = vtk.vtkFloatArray()
      vtkcoords.SetNumberOfComponents(3)
      vtkcoords.SetNumberOfTuples(self.nbDof)
      
      value = vtk.vtkDoubleArray()
      value.SetName("Value")
      
      for i in range(self.nbDof):
         vtkcoords.SetTuple3(i, self.meshNode[i,0], self.meshNode[i,1], self.meshNode[i,2])
         value.InsertNextValue(val[i])
      
      vtkpts = vtk.vtkPoints()
      vtkpts.SetData(vtkcoords)
      
      strips = vtk.vtkCellArray()
      for i in range(self.nbElem):
         strips.InsertNextCell(3)
         for j,inode in enumerate(self.meshTable[i,:]):
            strips.InsertCellPoint(inode)
            
      # Create the dataset. In this case, we create a vtkPolyData
      polydata = vtk.vtkPolyData()
      # Assign points and cells
      polydata.SetPoints(vtkpts)
      polydata.SetStrips(strips)
      polydata.GetPointData().SetScalars(value)
      
      return polydata

   def createVTKisocontour(self, dataVTK, isovalue, filename, poidsVTK, s):
      "Creation des isocontours"
   
      # Creation des contours
      contours = vtk.vtkContourFilter()
      contours.SetInput(dataVTK)
      #~ contours.GenerateValues(nb_isocontour, values.min(), values.max())
      contours.SetValue(0, isovalue)
      contours.Update()
      
      # Gestion des contours
      connectFilter = vtk.vtkPolyDataConnectivityFilter()
      connectFilter.SetInputConnection(contours.GetOutputPort())
      connectFilter.ScalarConnectivityOff() 
      connectFilter.ColorRegionsOn()
      connectFilter.SetExtractionModeToSpecifiedRegions()
      connectFilter.AddSpecifiedRegion(0)
      connectFilter.Update()
   
      if connectFilter.GetNumberOfExtractedRegions() == 1 :
         
         # Recuperation des ponderations
         probe = vtk.vtkProbeFilter()
         probe.SetInputConnection(connectFilter.GetOutputPort())
         probe.SetSource( poidsVTK)
         probe.Update()
         ponderation = VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray('Value'))
   
         numPoints = probe.GetOutput().GetNumberOfPoints() # get the number of points on the line
         #intialise the points on the line    
         x = np.zeros(numPoints)
         y = np.zeros(numPoints)
         z = np.zeros(numPoints)
         pts = np.zeros((numPoints , 3))
         #get the coordinates of the points on the line
         for i in range(numPoints):
            x[i],y[i],z[i] = probe.GetOutput().GetPoint(i)
         pts[:,0]=x
         pts[:,1]=y
         pts[:,2]=z
              
         # Organisation des noeuds
         pts,ponderation = self.orderingNodesOnLine(pts,ponderation)
         #~ plot3Dcontour(pts)
         
         # scaling ponderation
         ponderation, seuil = self.seuilMass(ponderation, s)
         
         # Estimation du centre de gravite
         center = self.getCdGfromContour(pts, ponderation, seuil)
   
         # Ecriture du contour pour post traitement
         #~ self.writeVTKContour(pts, filename)
         
      else :
         center = np.array([0.0,0.0,0.0,-1.0,0.0,0.0,0.0])
      
      return center,connectFilter.GetOutput()

   def orderingNodesOnLine(self, nodes, value):
      "Ordering nodes on lines"
   
      nb = nodes.shape[0]
      
      order_nodes = np.zeros(nodes.shape,float)
      order_value = np.zeros(value.shape,float)
      order_nodes[0,:] = nodes[0,:]
      order_value[0]   = value[0]
      nodes = np.delete(nodes,0,axis = 0)
      value = np.delete(value,0)
      
      
      for i in range(nb-2):
         dist = np.linalg.norm((nodes - order_nodes[i,:]),axis = 1)
         inds = dist.argsort()
         # Premier passage qui indique la direction de recherche
         #~ if i == 0 :
         order_nodes[1,:] = nodes[inds[0],:]
         order_value[1]   = value[inds[0]]
         vec = order_nodes[1,:] - order_nodes[0,:]
         nodes = np.delete(nodes,inds[0],axis = 0)
         value = np.delete(value,inds[0])
         # pour les autres segments on cherche le point le plus proches 
         # dans la bonne direction (produit scalaire > 0)
         #~ else :
            #~ ok = True
            #~ j = 0
            #~ while ok :
               #~ print j, inds[j], i
               #~ vec_j = nodes[inds[j],:] - order_nodes[i,:]
               #~ print np.dot(vec_j,vec)
               #~ ok = not(np.dot(vec_j,vec) > 0.0)
               #~ j += 1
            #~ #on a trouve le bon candidat
            #~ # le candidat j-1 est dans la bonne direction
            #~ order_nodes[i+1,:] = nodes[inds[j-1],:]
            #~ order_value[i+1]   = value[inds[j-1]]
            #~ vec = order_nodes[i+1,:] - order_nodes[i,:]
            #~ nodes = np.delete(nodes,inds[j-1],axis = 0)
            #~ value = np.delete(value,inds[j-1])
               
      # Il reste maintenant qu'un seul point que l'on ajoute
      order_nodes[nb-1,:] = nodes[0,:]
      order_value[nb-1]   = value[0]
      nodes = np.delete(nodes,0,axis = 0)
      value = np.delete(value,0)
         
      return order_nodes, order_value

   def getCdGfromLine(self, PtA, PtB, pdA, pdB, seuil):
      "get Inertie of 2 nodes lines type S2xxx"
      
      if (pdA > seuil) or (pdB > seuil) :
         AB = PtB - PtA
         L  = np.linalg.norm(AB)
         M  = L*((pdA + pdB)/2.0)
         
         OG = PtA + (((2.0/3.0)*pdB + (1.0/3.0)*pdA)*AB)/(pdA + pdB)
      else :
         print "Un segment n a pas de poids "
         OG = np.array([0.0, 0.0 ,0.0])
         M  = 0.0
      
      return OG, M

   def getCdGfromContour(self, pts, ponderation, seuil):
      "Calcul du centre de gravite d'un contour de segment "
      OG = np.zeros((3),float)
      M  = 0.0
      # Boucle sur les segments
      for i in range(pts.shape[0]-1):
         OGi, Mi = self.getCdGfromLine(pts[i,:], pts[i+1,:], ponderation[i], ponderation[i+1], seuil)
         OG = OG + Mi*OGi
         M  = M + Mi
      # Pour fermer le contour
      OGi, Mi = self.getCdGfromLine(pts[i+1,:], pts[0,:], ponderation[i+1], ponderation[0], seuil)
      OG = OG + Mi*OGi
      M  = M + Mi
      # Calcul du CdG et des proprietes geometriques
      CdG = OG/M
      R  = np.linalg.norm((pts - CdG),axis = 1)
      # Sauvegarde des proprietes
      center = np.concatenate((CdG,np.array([seuil, R.max(), R.min(), R.mean()])),axis = 1)
   
      return center


   def movingAverage(self, val, window_size):
      "Moyenne mobile"
      window= np.ones(int(window_size))/float(window_size)
      val[:,0] = np.convolve(val[:,0], window, 'same')
      val[:,1] = np.convolve(val[:,1], window, 'same')
      val[:,2] = np.convolve(val[:,2], window, 'same')
      return val[window_size:-window_size]
   
   def cutCenterLine(self,centerLine):
      "Supprime les parties non defini de la center line"
      ind = (centerLine[:,3]>=0.0).nonzero()[0]
      
      return centerLine[ind,:]

   def seuilMass(self, values, seuil):
      "Seuillage des data"
      s = values.mean()/seuil
   
      ind = (values < s).nonzero()[0]
      
      if len(ind)>0 :
         values[ind] = 0.0
   
      return values, s

   def removeBadContour(self, centerLine):
      "Supprime les mauvais contours"
      
      inds = (centerLine[:,3]<>-1.0).nonzero()[0]
      centerLine = centerLine[inds,:]
      return centerLine

   def writeVTKCenterLine(self, filename):
      "Create data in VTK format to compute isocontour"
      
      writer = vtk.vtkXMLPolyDataWriter()
      writer.SetFileName(filename+'.vtp')
      writer.SetInput(self.profileTubes.GetOutput())
      writer.Write()

   def createVTKContour(self, nodes):
      "charge un contour en memoire"

      # Conteneur du profil
      pts = vtk.vtkPoints()
      for n,node in enumerate(nodes):
         pts.InsertPoint(n,node[0],node[1],node[2])
      
      lines = vtk.vtkCellArray()
      lines.InsertNextCell(nodes.shape[0] + 1)
      for n,node in enumerate(nodes):
         lines.InsertCellPoint(n)
      lines.InsertCellPoint(0)
         
      self.profileContour = vtk.vtkPolyData()
      self.profileContour.SetPoints(pts)
      self.profileContour.SetLines(lines)
      
      # Add thickness to the resulting line.
      self.profileContourTubes = vtk.vtkTubeFilter()
      self.profileContourTubes.SetNumberOfSides(10)
      self.profileContourTubes.SetInput(self.profileContour)
      self.profileContourTubes.SetRadius(0.2)
      
   def writeVTKContour(self, filename):
      "Create data in VTK format to compute isocontour"
      # Add thickness to the resulting line.
      writer = vtk.vtkXMLPolyDataWriter()
      writer.SetFileName(filename+'.vtp')
      writer.SetInput(self.profileContour)
      writer.Write()

   def computeCenterLine(self):
      "Calcule la centerline"
      
      polydata0 = self.createVTKpolydata(self.X)
      polydata1 = self.createVTKpolydata(self.Y)
      s        = np.linspace(-0.99,0.99,self.resolution)
      self.centerLine = np.zeros((s.shape[0],7),float)
      for i,x in enumerate(s) :
         self.centerLine[i,:],poly = self.createVTKisocontour(polydata0, x, 'Contour_1_'+str(i), polydata1, 10.0)
         
      self.centerLine = self.removeBadContour(self.centerLine)
      self.centerLine = self.movingAverage(self.centerLine, self.filter)

   def createVTKCenterLine(self, type ):
      "Affiche la center line dans la fenetre graphique"
      
      nodes = self.centerLine
      
      # Conteneur du tube
      self.pts = vtk.vtkPoints()
      for n,node in enumerate(nodes):
         self.pts.InsertPoint(n,node[0],node[1],node[2])
      
      self.lines = vtk.vtkCellArray()
      self.lines.InsertNextCell(nodes.shape[0])
      for n,node in enumerate(nodes):
         self.lines.InsertCellPoint(n)
         
      self.scalars = vtk.vtkFloatArray()
      self.scalars.SetNumberOfValues(nodes.shape[0])
      for n,node in enumerate(nodes):
         if type == 'max' : self.scalars.SetValue(n, node[4])
         if type == 'min' : self.scalars.SetValue(n, node[5])
         if type == 'moy' : self.scalars.SetValue(n, node[6])
      
      self.profileData = vtk.vtkPolyData()
      self.profileData.SetPoints(self.pts)
      self.profileData.SetLines(self.lines)
      self.profileData.GetPointData().SetScalars(self.scalars)
      
      # Add thickness to the resulting line.
      self.profileTubes = vtk.vtkTubeFilter()
      self.profileTubes.SetNumberOfSides(10)
      self.profileTubes.SetInput(self.profileData)
      self.profileTubes.SetRadius(0.1)
      self.profileTubes.SetVaryRadiusToVaryRadiusByScalar()
      

      
      
      
      
      

