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
#! requirements in conditions enabling the security of their systems and/or 
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
import vtk
import numpy as np
import os
from vtk.util.numpy_support import vtk_to_numpy as vtk_to_np 
#~ import pylab as pl
#----------------------------------------------------------------------------------------------  
class vtkPipeline():
   """ Classe de fonction pour le pipe vtk """
   
   def __init__(self, widget, ren, config):
      """ Initialisation des variables locales de la classe """
      self.widget = widget
      self.ren    = ren
      
      self.path = []
      self.listFiles = []
      self.col = 0
      self.lig = 0
      self.nb_coupe = 0
      self.SizeOfVolum = np.array([0.0, 0.0, 0.0])
      self.pixel_spacing_x = 0.0
      self.pixel_spacing_y = 0.0
      self.slice_spacing = 0.0
      self.min_pixel_value = 0
      self.max_pixel_value = 0
      self.GeoInfos = ""
      self.PatientInfos = ""
      self.valTransversal = 0
      self.valCoronal = 0
      self.valSagittal = 0
      self.valSeuil = int(config['VALSEUIL'])
      self.valFilter = int(config['VALFILTER'])
      self.valGaussianFilter = int(config['VALGAUSSIANFILTER'])
      self.valOpacity = 10
      self.valSmoothIter = int(config['VALSMOOTHITER'])
      self.valDecimate = float(config['VALDECIMATE'])/100.0
      self.VOI_pt1 = [0.0,0.0,0.0]
      self.VOI_pt2 = [0.0,0.0,0.0]
      self.SurfaceArea = 0.0
      self.SurfaceVolum = 0.0
      
      self.aneurysmExtractor = None
      self.aneurysm = None
      self.outlineVOI = None
      self.smoothAneurysm = None
      self.sagittalCut = None
      self.coronalCut = None
      self.transversalCut = None
      self.Axis = None
      self.outline = None
      self.vtkInfoGeoTitle = None
      self.vtkInfoGeo = None
      self.vtkInfoPatTitle = None
      self.vtkInfoPat = None
      self.meshActor = None
      self.mesh3DActor = None
      self.outputMesh = None
      self.output3DMesh = None
      self.contourActor = None
      self.centerlineActor = None


   def vtkCreateGeometricInformations(self):
      """ Creation du texte d'information geometrique"""
      self.GeoInfos = "--------------------------------\n"
      self.GeoInfos += "nb pixel x : %i\n" %self.lig
      self.GeoInfos += "nb pixel y : %i\n" %self.col
      self.GeoInfos += "nb slice z : %i\n" %self.nb_coupe
      self.GeoInfos += "--------------------------------\n"
      self.GeoInfos += "spacing x : %5.3f\n" %self.pixel_spacing_x
      self.GeoInfos += "spacing y : %5.3f\n" %self.pixel_spacing_y
      self.GeoInfos += "spacing z : %5.3f\n" %self.slice_spacing
      self.GeoInfos += "--------------------------------\n"
      self.GeoInfos += "length x : %5.3f\n" %self.SizeOfVolum[0]
      self.GeoInfos += "length y : %5.3f\n" %self.SizeOfVolum[1]
      self.GeoInfos += "length z : %5.3f\n" %self.SizeOfVolum[2]
      self.GeoInfos += "--------------------------------\n"
      self.GeoInfos += "minimal intensity : %i\n" %self.min_pixel_value 
      self.GeoInfos += "maximal intensity : %i\n" %self.max_pixel_value 
      self.GeoInfos += "--------------------------------\n"
      self.GeoInfos += "surface intensity = %i\n" %self.valSeuil
      self.GeoInfos += "Surface area  = %5.3f\n" %self.SurfaceArea
      self.GeoInfos += "Surface volum = %5.3f\n" %self.SurfaceVolum
      self.GeoInfos += "--------------------------------\n"

   def vtkInitialData(self, read_path, volum = 'GAUSSIAN'):
      """ Creation des donnees initiales """
      self.path = read_path
      #-----------------------------------------
      # Creation de la donnee volumique
      #-----------------------------------------
      self.vtkReadVolumAll(os.path.join(self.path,'image.raw'))
      self.valSagittal = (self.lig-1)/2
      self.valCoronal = (self.col-1)/2
      self.valTransversal = (self.nb_coupe-1)/2
      self.min_pixel_value = self.NumpyData.min()
      self.max_pixel_value = self.NumpyData.max()
      if self.valSeuil >= self.max_pixel_value :
         self.valSeuil = (self.max_pixel_value - self.min_pixel_value)/2
      self.VOI_pt2 = [self.col,self.lig,self.nb_coupe]
      #-----------------------------------------
      # Start by creatin a black/white lookup table.
      #-----------------------------------------
      self.bwLut = vtk.vtkLookupTable()
      self.bwLut.SetTableRange(self.min_pixel_value, self.max_pixel_value)
      self.bwLut.SetSaturationRange(0, 0)
      self.bwLut.SetHueRange(0, 0)
      self.bwLut.SetValueRange(0, 1)
      self.bwLut.Build()
      #-----------------------------------------
      # Creation de la bounding box
      #-----------------------------------------
      self.vtkBoundingBox()
      #-----------------------------------------
      # Creation de la VOI
      #-----------------------------------------
      self.vtkVOIdata(volum = volum)
      self.vtkBoundingVOI()
      #-----------------------------------------
      # Creation des infos geometriques
      #-----------------------------------------
      self.vtkCreateGeometricInformations()
      #-----------------------------------------
      # Creation de la couleur des coupe
      #-----------------------------------------
      self.sliceColors0 = vtk.vtkImageMapToColors()
      self.sliceColors0.SetInputConnection(self.vtkVolumAll.GetOutputPort())
      self.sliceColors0.SetLookupTable(self.bwLut)
      self.sliceColors1 = vtk.vtkImageMapToColors()
      self.sliceColors1.SetInputConnection(self.vtkVolumAll.GetOutputPort())
      self.sliceColors1.SetLookupTable(self.bwLut)
      self.sliceColors2 = vtk.vtkImageMapToColors()
      self.sliceColors2.SetInputConnection(self.vtkVolumAll.GetOutputPort())
      self.sliceColors2.SetLookupTable(self.bwLut)
      #-----------------------------------------
      # Creation de la camera
      #-----------------------------------------
      self.vtkCreateCamera()
      #-----------------------------------------
      # Fin de creation des elements et affichage
      #-----------------------------------------
      self.vtkRedraw()

   def vtkComputeGeometricSurface(self):
      """ Compute th area and volum of surface construction"""
      measured_polydata = vtk.vtkMassProperties()
      measured_polydata.SetInput(self.smoothAneurysm.GetOutput())
      self.SurfaceArea = measured_polydata.GetSurfaceArea()
      self.SurfaceVolum = measured_polydata.GetVolume()

   def vtkMoveVOI(self):
      """ Extraction du volume d'interet"""
      #-----------------------------------------
      # Mouvement du VOI
      #-----------------------------------------
      self.vtkVolum.SetVOI(int(self.VOI_pt1[0]), int(self.VOI_pt2[0])  , \
                           int(self.VOI_pt1[1]), int(self.VOI_pt2[1])  , \
                           int(self.VOI_pt1[2]), int(self.VOI_pt2[2]))
      self.vtkRedraw()
   
   def vtkCreateIsoContour(self , config = 'MARCHINGCUBES'):
      """ Fonction pour la creation de l'isocontour pour une valeur de seuillage"""
      
      #-----------------------------------------
      # Creation de l isocontour
      #-----------------------------------------
      if config == 'CONTOUR' :
         self.aneurysmExtractor = vtk.vtkContourFilter()
      if config == 'MARCHINGCUBES' :
         self.aneurysmExtractor = vtk.vtkMarchingCubes()
         
      self.aneurysmExtractor.SetInputConnection(self.vtkVolumBlur.GetOutputPort())
      self.aneurysmExtractor.SetValue(0, self.valSeuil)
      
      self.aneurysmExtractor.ComputeNormalsOn()
         
      self.aneurysmTriangleFilter = vtk.vtkTriangleFilter()
      self.aneurysmTriangleFilter.SetInputConnection(self.aneurysmExtractor.GetOutputPort())
  
      self.aneurysmCleanFilter = vtk.vtkCleanPolyData()
      self.aneurysmCleanFilter.SetInputConnection(self.aneurysmTriangleFilter.GetOutputPort()) 
      
      self.aneurysmConnectFilter = vtk.vtkPolyDataConnectivityFilter()
      self.aneurysmConnectFilter.SetExtractionModeToLargestRegion()
      self.aneurysmConnectFilter.ScalarConnectivityOff() 
      self.aneurysmConnectFilter.SetInputConnection(self.aneurysmCleanFilter.GetOutputPort()) 
      
      self.aneurysmNormals = vtk.vtkPolyDataNormals()
      self.aneurysmNormals.SetInputConnection(self.aneurysmConnectFilter.GetOutputPort())
      self.aneurysmNormals.SetFeatureAngle(60.0)
      self.aneurysmNormals.ComputeCellNormalsOn()
      self.aneurysmNormals.ComputePointNormalsOn()
      self.aneurysmNormals.ConsistencyOn()
      self.aneurysmNormals.AutoOrientNormalsOn()
      self.aneurysmNormals.Update()
      
   def vtkDeleteIsoContour(self):
      """ Fonction de suppression du calcul du contour """
      self.aneurysmExtractor = None
      self.aneurysmTriangleFilter = None
      self.aneurysmCleanFilter = None
      self.aneurysmConnectFilter = None
      self.aneurysmNormals = None

   def vtkCreateFilterIsoContour(self):
      """ Fonction pour le filtrage du contour"""
      #-----------------------------------------
      # Application d'un filtre de subdivision
      #-----------------------------------------
      self.subdiviseAneurysm = vtk.vtkLoopSubdivisionFilter()
      self.subdiviseAneurysm.SetNumberOfSubdivisions(self.valFilter)
      self.subdiviseAneurysm.SetInputConnection(self.aneurysmNormals.GetOutputPort())
      
      #-----------------------------------------
      # Application d'un filtre decimate et smooth
      #-----------------------------------------
      self.aneurysmDecimator = vtk.vtkDecimatePro()
      self.aneurysmDecimator.SetInputConnection(self.subdiviseAneurysm.GetOutputPort()) 
      self.aneurysmDecimator.SetTargetReduction(self.valDecimate) 
      self.aneurysmDecimator.PreserveTopologyOn() 
      
      self.smoothAneurysm = vtk.vtkSmoothPolyDataFilter()
      self.smoothAneurysm.SetInputConnection(self.aneurysmDecimator.GetOutputPort()) 
      self.smoothAneurysm.SetNumberOfIterations(self.valSmoothIter) 
      self.smoothAneurysm.SetRelaxationFactor(0.07) 
      self.smoothAneurysm.SetFeatureAngle(90) 
      self.smoothAneurysm.FeatureEdgeSmoothingOn()
      #-----------------------------------------
      # Sauvegarde des datas pour le calcul des centerlines
      #-----------------------------------------
      surface_filter = vtk.vtkDataSetSurfaceFilter()
      surface_filter.SetInputConnection(self.smoothAneurysm.GetOutputPort())
      
      meshNormals = vtk.vtkPolyDataNormals()
      meshNormals.SetInputConnection(surface_filter.GetOutputPort())
      meshNormals.Update()
      
      self.surfNode  = vtk_to_np(self.smoothAneurysm.GetOutput().GetPoints().GetData())
      self.surfTable = vtk_to_np(self.smoothAneurysm.GetOutput().GetPolys().GetData())
      
      #-----------------------------------------
      # Realisation du mapper ppour affichage
      #-----------------------------------------
      self.aneurysmStripper = vtk.vtkStripper()
      self.aneurysmStripper.SetInputConnection(self.smoothAneurysm.GetOutputPort())
      
      self.aneurysmMapper = vtk.vtkPolyDataMapper()
      self.aneurysmMapper.SetInputConnection(self.aneurysmStripper.GetOutputPort())
      self.aneurysmMapper.ScalarVisibilityOff()
      self.aneurysmMapper.SetScalarRange(self.min_pixel_value, self.max_pixel_value)
      self.aneurysmMapper.SetScalarModeToUsePointData()
      self.aneurysmMapper.SetLookupTable(self.bwLut)
      
   def vtkDeleteBlockFilterIsoContour(self):
      """ Fonction de suppression du filtrage du contour"""
      self.subdiviseAneurysm = None
      self.aneurysmDecimator = None
      self.smoothAneurysm = None
      self.aneurysmStripper = None
      
   def vtkPlotIsoContour(self):
      """ Fonction de creation 3D de l'iso contour"""
      self.aneurysm = vtk.vtkActor()
      self.aneurysm.SetMapper(self.aneurysmMapper)
      self.aneurysm.GetProperty().SetDiffuseColor(1, 1, .9412)
      self.aneurysm.GetProperty().SetOpacity(self.valOpacity/10.0)
      self.ren.AddActor(self.aneurysm)
      self.vtkRedraw()
      
   def vtkDeletePlotIsoContour(self):
      """ Fonction de suppression de l'affichage"""
      self.ren.RemoveActor(self.aneurysm)
      self.vtkRedraw()
      self.aneurysm = None
      
   def vtkDeleteAllActor(self):
      """ Fonction de suppression des actors vtk """
      if self.aneurysm        <> None : self.ren.RemoveActor(self.aneurysm)
      if self.outlineVOI      <> None : self.ren.RemoveActor(self.outlineVOI)
      if self.smoothAneurysm  <> None : self.ren.RemoveActor(self.smoothAneurysm)
      if self.sagittalCut     <> None : self.ren.RemoveActor(self.sagittalCut)
      if self.coronalCut      <> None : self.ren.RemoveActor(self.coronalCut)
      if self.transversalCut  <> None : self.ren.RemoveActor(self.transversalCut)
      if self.Axis            <> None : self.ren.RemoveActor(self.Axis)
      if self.outline         <> None : self.ren.RemoveActor(self.outline)
      if self.vtkInfoGeoTitle <> None : self.ren.RemoveActor(self.vtkInfoGeoTitle)
      if self.vtkInfoGeo      <> None : self.ren.RemoveActor(self.vtkInfoGeo)
      if self.vtkInfoPatTitle <> None : self.ren.RemoveActor(self.vtkInfoPatTitle)
      if self.vtkInfoPat      <> None : self.ren.RemoveActor(self.vtkInfoPat)
      if self.meshActor       <> None : self.ren.RemoveActor(self.meshActor)
      self.vtkRedraw()
      self.aneurysm = None
      self.outlineVOI = None
      self.smoothAneurysm = None
      self.sagittalCut = None
      self.coronalCut = None
      self.transversalCut = None
      self.Axis = None
      self.outline = None
      self.vtkInfoGeoTitle = None
      self.vtkInfoGeo = None
      self.vtkInfoPatTitle = None
      self.vtkInfoPat = None
      self.meshActor = None

   def vtkReadVolumAll(self, filename):
      """ Chargement des fichiers images dicom"""
      #-----------------------------------------
      # Creation de la donnee volumique
      #-----------------------------------------
      self.vtkVolumAll = vtk.vtkVolume16Reader()
      self.vtkVolumAll.SetDataDimensions(self.col, self.lig)
      self.vtkVolumAll.SetDataByteOrderToLittleEndian()
      self.vtkVolumAll.SetFilePrefix(filename)
      self.vtkVolumAll.SetImageRange(1, self.nb_coupe)
      self.vtkVolumAll.SetDataSpacing(self.pixel_spacing_x, self.pixel_spacing_y, self.slice_spacing)
      self.SizeOfVolum = np.array([self.lig*self.pixel_spacing_x,self.col*self.pixel_spacing_y,self.slice_spacing*self.nb_coupe])
      
      self.NumpyData = np.zeros((self.lig,self.col,self.nb_coupe),int)
      for i in range(self.nb_coupe):
         self.NumpyData[:,:,i] = vtk_to_np(self.vtkVolumAll.GetImage(i+1).GetPointData().GetArray(0)).reshape(self.lig,self.col)
         
   def vtkVOIdata(self, volum = 'GAUSSIAN'):
      """ Modification du volume d'interet """
      #-----------------------------------------
      # Creation de la VOI
      #-----------------------------------------
      self.vtkVolum = vtk.vtkExtractVOI()
      self.vtkVolum.SetInput(self.vtkVolumAll.GetOutput())
      self.vtkVolum.SetVOI(self.VOI_pt1[0], self.VOI_pt2[0], \
                           self.VOI_pt1[1], self.VOI_pt2[1], \
                           self.VOI_pt1[2], self.VOI_pt2[2])
      
      self.vtkVolumBlur = vtk.vtkImageGaussianSmooth()
      #-----------------------------------------
      # Creation du volume sur volume filtre gaussien
      #-----------------------------------------
      if volum == 'GAUSSIAN':
         self.vtkVolumBlur.SetInputConnection(self.vtkVolum.GetOutputPort()) 
      #-----------------------------------------
      # Creation du volume sur volume filtre gradient + gaussien
      #-----------------------------------------
      if volum == 'GRADIENT':
         self.vtkVolumGrad = vtk.vtkImageGradientMagnitude()
         self.vtkVolumGrad.SetDimensionality(3)
         self.vtkVolumGrad.HandleBoundariesOn()
         self.vtkVolumGrad.SetInputConnection(self.vtkVolum.GetOutputPort())
         self.vtkVolumBlur.SetInputConnection(self.vtkVolumGrad.GetOutputPort()) 
      
      self.vtkVolumBlur.SetStandardDeviation(self.valGaussianFilter)
      self.vtkVolumBlur.Update()
      self.vtkVolumRange = self.vtkVolumBlur.GetOutput().GetScalarRange()
      
      print " Range from Volume Of Interest : ", self.vtkVolumRange
      
   def vtkHistogram(self):
      """Fonction de creation de l'histogramme """
      #-----------------------------------------
      # Creation de l'histogramme
      #-----------------------------------------
      data = self.NumpyData[int(self.VOI_pt1[1]):int(self.VOI_pt2[1]),\
                            int(self.VOI_pt1[0]):int(self.VOI_pt2[0]),\
                            int(self.VOI_pt1[2]):int(self.VOI_pt2[2])]
      
      data = data.reshape(data.shape[0]*data.shape[1]*data.shape[2],1)
      
      print "Maximum : ",data.max()
      print "Minimum : ",data.min()
      
      Histogram = np.histogram(data,np.arange(data.min(),data.max(),1))

   def vtkCreateCamera(self):
      """ Creation de la camera de visualisation """
      #-----------------------------------------
      # Creation de la camera
      #-----------------------------------------
      self.vtkCam = vtk.vtkCamera()
      self.vtkCam.SetPosition(self.SizeOfVolum[0]/2.0, self.SizeOfVolum[1]/2.0, self.SizeOfVolum.max()*4)
      self.vtkCam.SetFocalPoint(self.SizeOfVolum[0]/2.0, self.SizeOfVolum[1]/2.0, 0.0)
      self.vtkCam.SetViewUp(0, 1, 0)
      self.ren.SetActiveCamera(self.vtkCam)
      self.ren.ResetCamera()
   
   def vtkCameraViewZ(self):
      """ Place la camera de normale Z """
      self.vtkCam.SetPosition(self.SizeOfVolum[0]/2.0, self.SizeOfVolum[1]/2.0, self.SizeOfVolum.max()*4)
      self.vtkCam.SetFocalPoint(self.SizeOfVolum[0]/2.0, self.SizeOfVolum[1]/2.0, 0.0)
      self.vtkCam.SetViewUp(0, 1, 0)
      self.vtkRedraw()
   
   def vtkCameraViewY(self):
      """ Place la camera de normale Y """
      self.vtkCam.SetPosition(self.SizeOfVolum[0]/2.0, self.SizeOfVolum.max()*4, self.SizeOfVolum[2]/2.0)
      self.vtkCam.SetFocalPoint(self.SizeOfVolum[0]/2.0, 0.0, self.SizeOfVolum[2]/2.0)
      self.vtkCam.SetViewUp(0, 0, 1)
      self.vtkRedraw()
      
   def vtkCameraViewX(self):
      """ Place la camera de normale X """
      self.vtkCam.SetPosition(self.SizeOfVolum.max()*4, self.SizeOfVolum[1]/2.0, self.SizeOfVolum[2]/2.0)
      self.vtkCam.SetFocalPoint(0.0, self.SizeOfVolum[1]/2.0, self.SizeOfVolum[2]/2.0)
      self.vtkCam.SetViewUp(0, 1, 0)
      self.vtkRedraw()
   
   def vtkAxis(self):
      """ Creation des axes"""
      transform = vtk.vtkTransform()
      transform.Scale(0.1*self.SizeOfVolum.max(), 0.1*self.SizeOfVolum.max(), 0.1*self.SizeOfVolum.max())
      self.Axis = vtk.vtkAxesActor()
      self.Axis.SetUserTransform(transform)
      self.ren.AddActor(self.Axis)
      self.vtkRedraw()
   
   def vtkDeleteAxis(self):
      """ Effacement des axes """
      self.ren.RemoveActor(self.Axis)
      self.vtkRedraw()
      self.Axis = None
   
   def vtkBoundingBox(self):
      """ Creation du  polydata de la bounding box """
      self.outlineData = vtk.vtkOutlineFilter()
      self.outlineData.SetInputConnection(self.vtkVolumAll.GetOutputPort())
      self.mapOutline = vtk.vtkPolyDataMapper()
      self.mapOutline.SetInputConnection(self.outlineData.GetOutputPort())

   def vtkPlotBoundingBox(self):
      """ Creation visuelle de la bounding box"""
      self.outline = vtk.vtkActor()
      self.outline.SetMapper(self.mapOutline)
      self.outline.GetProperty().SetColor(1, 1, 1)
      self.ren.AddActor(self.outline)
      self.vtkRedraw()
      
   def vtkDeleteBoundingBox(self):
      """ Effacement de la bounding box """
      self.ren.RemoveActor(self.outline)
      self.vtkRedraw()
      self.outline = None

   def vtkBoundingVOI(self):
      """ Creation du  polydata de la bounding VOI"""
      self.outlineDataVOI = vtk.vtkOutlineFilter()
      self.outlineDataVOI.SetInputConnection(self.vtkVolum.GetOutputPort())
      self.mapOutlineVOI = vtk.vtkPolyDataMapper()
      self.mapOutlineVOI.SetInputConnection(self.outlineDataVOI.GetOutputPort())
   
   def vtkPlotBoundingVOI(self):
      """ Affichage de la bounding VOI"""
      self.outlineVOI = vtk.vtkActor()
      self.outlineVOI.SetMapper(self.mapOutlineVOI)
      self.outlineVOI.GetProperty().SetColor(0.95, 0.18, 0.18)
      self.ren.AddActor(self.outlineVOI)
      self.vtkRedraw()
      
   def vtkDeleteBoundingVOI(self):
      """ Effacement de la bounding box """
      self.ren.RemoveActor(self.outlineVOI)
      self.vtkRedraw()
      self.outlineVOI = None

   def vtkCreateSTLContour(self):
      """ Affichage d'un fichier STL"""
      self.STLaneurysmMapper = vtk.vtkPolyDataMapper()
      self.smoothAneurysm = None
      self.smoothAneurysm = vtk.vtkActor()
      self.smoothAneurysm.SetMapper(self.STLaneurysmMapper)
      self.smoothAneurysm.GetProperty().SetDiffuseColor(1, 1, .9412)
      self.smoothAneurysm.GetProperty().SetOpacity(self.valOpacity/10.0)
      self.ren.AddActor(self.smoothAneurysm)
      self.vtkRedraw()
      
   def vtkCreateSTLMesh(self):
      """ Affichage d'un fichier mesh STL"""
      
      self.meshMapper = vtk.vtkReverseSense()
      #~ self.meshMapper.SetInputConnection(self.meshNormals.GetOutputPort())
      self.meshMapper.ReverseCellsOn()
      self.meshMapper.ReverseNormalsOn()
      
      self.STLMeshMapper = vtk.vtkPolyDataMapper()
      self.STLMeshMapper.SetInputConnection(self.meshMapper.GetOutputPort())
      
      self.meshNormals = vtk.vtkPolyDataNormals()
      self.meshNormals.SetInputConnection(self.STLMeshMapper.GetOutputPort())
      self.meshNormals.SetFeatureAngle(60.0)
      self.meshNormals.ComputeCellNormalsOn()
      self.meshNormals.ComputePointNormalsOn()
      self.meshNormals.ConsistencyOn()
      self.meshNormals.AutoOrientNormalsOn()

   def vtkPlotSTLMesh(self):
      "Creation du plot STLMesh"
      self.meshActor = None
      self.meshActor = vtk.vtkActor()
      self.meshActor.SetMapper(self.meshMapper)
      #~ self.meshActor.GetProperty().SetDiffuseColor(1, 0.2, 0.2)
      #~ self.meshActor.GetProperty().SetOpacity(1.0)
      
      self.meshActor.GetProperty().SetDiffuseColor(1.0000, 0.3882, 0.2784)
      self.meshActor.GetProperty().SetSpecularColor(1, 1, 1)
      self.meshActor.GetProperty().SetSpecular(.4)
      self.meshActor.GetProperty().SetSpecularPower(50)
      
      self.ren.AddActor(self.meshActor)
      self.vtkRedraw()
      
   def vtkImportSTLfile(self, path):
      """ Importation d'un fichier stl """
      stlReader = vtk.vtkSTLReader()
      stlReader.SetFileName(path)
      stlReader.Update()
      
      reverse = vtk.vtkReverseSense()
      reverse.SetInputConnection(stlReader.GetOutputPort())
      reverse.ReverseCellsOn()
      reverse.ReverseNormalsOn()
      
      self.meshMapper = vtk.vtkPolyDataMapper()
      if vtk.VTK_MAJOR_VERSION <= 5:
          self.meshMapper.SetInput(reverse.GetOutput())
      else:
          self.meshMapper.SetInputConnection(reverse.GetOutputPort())
      self.vtkPlotSTLMesh()

   def vtkWriteSTLfile(self, path, type, entity_type):
      """ Ecriture d'un fichier stl"""
      stlWriter = vtk.vtkSTLWriter()
      if type == 'ASCII'  : stlWriter.SetFileTypeToASCII()
      if type == 'BINARY' : stlWriter.SetFileTypeToBinary()
      stlWriter.SetFileName(path)
      if entity_type == 'Surface' : 
         stlWriter.SetInputConnection(self.smoothAneurysm.GetOutputPort())
      if entity_type == 'Mesh'    :
         surface_filter = vtk.vtkDataSetSurfaceFilter()
         surface_filter.SetInput(self.outputMesh)
         clean = vtk.vtkCleanPolyData()
         clean.SetInputConnection(surface_filter.GetOutputPort())
         clean.SetTolerance(0.001)
         triangle_filter = vtk.vtkTriangleFilter()
         triangle_filter.SetInputConnection(surface_filter.GetOutputPort())
         stlWriter.SetInputConnection(triangle_filter.GetOutputPort())
      stlWriter.Write()
      
   def vtkReadVTKfile(self, filename):
      """ Lecture d'un fichier vtk"""
      # Read the source file.
      reader = vtk.vtkUnstructuredGridReader()
      reader.SetFileName(filename)
      reader.Update() # Needed because of GetScalarRange
      
      surface_filter = vtk.vtkDataSetSurfaceFilter()
      surface_filter.SetInputConnection(reader.GetOutputPort()) 
      
      self.meshNormals = vtk.vtkPolyDataNormals()
      self.meshNormals.SetInputConnection(surface_filter.GetOutputPort())
      self.meshNormals.SetFeatureAngle(60.0)
      self.meshNormals.ComputeCellNormalsOn()
      self.meshNormals.ComputePointNormalsOn()
      self.meshNormals.ConsistencyOn()
      self.meshNormals.AutoOrientNormalsOn()
      self.meshNormals.Update()
      
      # Sauvegarde des proprietes maillage GMSH pour calcul centerline
      self.gmeshNode  = vtk_to_np(surface_filter.GetOutput().GetPoints().GetData())
      self.gmeshTable = vtk_to_np(surface_filter.GetOutput().GetPolys().GetData())
      
      # The quadric has nasty discontinuities from the way the edges are generated
      # so let's pass it though a CleanPolyDataFilter and merge any points which
      # are coincident, or very close
 
      self.outputMesh = vtk.vtkCleanPolyData()
      self.outputMesh.SetInputConnection(surface_filter.GetOutputPort())
      self.outputMesh.SetTolerance(0.001)
      
   def vtkPlotMesh(self, courbure = "MEAN"):
      """ Creation du maillage surfacique """
      self.meshActor = None
      self.vtkMeshCurvature = None
      # Creation de l'information des courbures
      self.vtkMeshCurvature = vtk.vtkCurvatures()
      if courbure == "GAUSS" : self.vtkMeshCurvature.SetCurvatureTypeToGaussian()
      if courbure == "MEAN"  : self.vtkMeshCurvature.SetCurvatureTypeToMean()
      self.vtkMeshCurvature.SetInputConnection(self.outputMesh.GetOutputPort())
      self.vtkMeshCurvature.Update()
      range_curvature = self.vtkMeshCurvature.GetOutput().GetScalarRange()

      # Create color of curvature
      Lut = vtk.vtkLookupTable()
      Lut.SetNumberOfColors(256)
      Lut.SetRange(range_curvature[0], range_curvature[1])
      Lut.Build()
      # Create the mapper that corresponds the objects of the vtk file
      # into graphics elements
      self.meshMapper = vtk.vtkPolyDataMapper()
      self.meshMapper.SetInputConnection(self.meshNormals.GetOutputPort())
      self.meshMapper.SetLookupTable(Lut)
      self.meshMapper.SetUseLookupTableScalarRange(1)
      # Create the Actor
      self.vtkDeleteMesh()
      self.meshActor = vtk.vtkActor()
      self.meshActor.SetMapper(self.meshMapper)
      self.meshActor.GetProperty().SetDiffuseColor(1, 0.2, 0.2)
      self.ren.AddActor(self.meshActor)
      self.vtkRedraw()
      
   def vtkDeleteMesh(self):
      """ Suppression du maillage surfacique """
      self.ren.RemoveActor(self.meshActor)
      self.vtkRedraw()
      self.meshActor = None
      self.vtkMeshCurvature = None

   def vtkReadVTKSolidfile(self, filename):
      """ Lecture d'un fichier vtk"""
      # Read the source file.
      reader = vtk.vtkUnstructuredGridReader()
      reader.SetFileName(filename)
      reader.Update() # Needed because of GetScalarRange
      output = reader.GetOutput()
      # Create the mapper that corresponds the objects of the vtk file
      # into graphics elements
      self.mesh3DMapper = vtk.vtkDataSetMapper()
      self.mesh3DMapper.SetInput(output)
      # Create the Actor
      self.mesh3DActor = vtk.vtkActor()
      self.mesh3DActor.SetMapper(self.mesh3DMapper)
      self.mesh3DActor.GetProperty().SetDiffuseColor(0.0, 0.8, 0.0)
      self.ren.AddActor(self.mesh3DActor)
      self.vtkRedraw()

   def vtkDeleteSolidMesh(self):
      """ Suppression du maillage 3D """
      self.ren.RemoveActor(self.mesh3DActor)
      self.vtkRedraw()
      self.mesh3DActor = None

   def vtkInfoPatient(self):
      """ Creation des infos patient """
      #-----------------------------------------
      # Creation de texte vtk
      #-----------------------------------------
      self.vtkInfoPatTitle = vtk.vtkTextActor()
      self.vtkInfoPatTitle.SetInput("Patient informations")
      txtprop=self.vtkInfoPatTitle.GetTextProperty()
      txtprop.SetFontFamilyToArial()
      txtprop.SetFontSize(14)
      txtprop.SetColor(1,1,1)
      self.vtkInfoPatTitle.SetDisplayPosition(10,180)
      #-----------------------------------------
      self.vtkInfoPat = vtk.vtkTextActor()
      self.vtkInfoPat.SetInput(self.PatientInfos)
      txtprop=self.vtkInfoPat.GetTextProperty()
      txtprop.SetFontFamilyToArial()
      txtprop.SetFontSize(12)
      txtprop.SetColor(1,1,1)
      self.vtkInfoPat.SetDisplayPosition(10,10)
      #-----------------------------------------
      self.ren.AddActor(self.vtkInfoPatTitle)
      self.ren.AddActor(self.vtkInfoPat)
      self.vtkRedraw()
   
   def vtkDeleteInfoPatient(self):
      """ Suppression des infos patient """
      self.ren.RemoveActor(self.vtkInfoPat)
      self.ren.RemoveActor(self.vtkInfoPatTitle)
      self.vtkRedraw()
      self.vtkInfoPat = None
      self.vtkInfoPatTitle = None
      
   def vtkInfoGeometric(self):
      """ Creation des infos geometriques """
      #-----------------------------------------
      # Creation de texte vtk
      #-----------------------------------------
      self.vtkInfoGeoTitle = vtk.vtkTextActor()
      self.vtkInfoGeoTitle.SetInput("Geometric informations :")
      txtprop=self.vtkInfoGeoTitle.GetTextProperty()
      txtprop.SetFontFamilyToArial()
      txtprop.SetFontSize(14)
      txtprop.SetColor(1,1,1)
      self.vtkInfoGeoTitle.SetDisplayPosition(10,280)
      #-----------------------------------------
      self.vtkInfoGeo = vtk.vtkTextActor()
      self.vtkInfoGeo.SetInput(self.GeoInfos)
      txtprop=self.vtkInfoGeo.GetTextProperty()
      txtprop.SetFontFamilyToArial()
      txtprop.SetFontSize(12)
      txtprop.SetColor(1,1,1)
      self.vtkInfoGeo.SetDisplayPosition(10,10)
      #-----------------------------------------
      self.ren.AddActor(self.vtkInfoGeoTitle)
      self.ren.AddActor(self.vtkInfoGeo)
      self.vtkRedraw()
   
   def vtkDeleteInfoGeometric(self):
      """ Suppression des infos geo"""
      self.ren.RemoveActor(self.vtkInfoGeo)
      self.ren.RemoveActor(self.vtkInfoGeoTitle)
      self.vtkRedraw()
      self.vtkInfoGeo = None
      self.vtkInfoGeoTitle = None
   
   def vtkTransversalSlice(self):
      """ Creation de la coupe transversale """
      self.transversalCut = vtk.vtkImageActor()
      self.transversalCut.SetInput(self.sliceColors2.GetOutput())
      self.transversalCut.SetDisplayExtent(0, self.col-1, 0, self.lig-1, self.valTransversal, self.valTransversal)
      #-----------------------------------------
      self.ren.AddActor(self.transversalCut)
      self.vtkRedraw()
      
   def vtkDeleteTransversalSlice(self):
      """ Suppression de la coupe transversal """
      self.ren.RemoveActor(self.transversalCut)
      self.vtkRedraw()
      self.transversalCut = None
      
   def vtkCoronalSlice(self):
      """ Creation de la coupe coronale """
      self.coronalCut = vtk.vtkImageActor()
      self.coronalCut.SetInput(self.sliceColors1.GetOutput())
      self.coronalCut.SetDisplayExtent(self.valCoronal, self.valCoronal, 0, self.lig-1, 0, self.nb_coupe - 1)
      #-----------------------------------------
      self.ren.AddActor(self.coronalCut)
      self.vtkRedraw()

   def vtkDeleteCoronalSlice(self):
      """ Suppression de la coupe transversal """
      self.ren.RemoveActor(self.coronalCut)
      self.vtkRedraw()
      self.coronalCut = None

   def vtkSagittalSlice(self):
      """ Creation de la coupe coronale """
      self.sagittalCut = vtk.vtkImageActor()
      self.sagittalCut.SetInput(self.sliceColors0.GetOutput())
      self.sagittalCut.SetDisplayExtent(0, self.col-1, self.valSagittal, self.valSagittal, 0, self.nb_coupe - 1)
      #-----------------------------------------
      self.ren.AddActor(self.sagittalCut)
      self.vtkRedraw()

   def vtkDeleteSagittalSlice(self):
      """ Suppression de la coupe transversal """
      self.ren.RemoveActor(self.sagittalCut)
      self.vtkRedraw()
      self.sagittalCut = None

   def vtkRedraw(self):
      """ re Draw the scene. """
      self.widget.Render()
      
   def vtkPlotCenterline(self, vtkCenterLine):
      "Affichage de la centerline"
      
      self.ren.RemoveActor(self.centerlineActor)
      self.vtkRedraw()
      self.centerlineActor = None
      
      reverse = vtk.vtkReverseSense()
      reverse.SetInputConnection(vtkCenterLine.GetOutputPort())
      reverse.ReverseCellsOff()
      reverse.ReverseNormalsOff()
      
      self.centerLineMapper = vtk.vtkPolyDataMapper()
      if vtk.VTK_MAJOR_VERSION <= 5:
          self.centerLineMapper.SetInput(reverse.GetOutput())
      else:
          self.centerLineMapper.SetInputConnection(reverse.GetOutputPort())
      
      self.centerlineActor = vtk.vtkActor()
      self.centerlineActor.SetMapper(self.centerLineMapper)
      self.centerlineActor.GetProperty().SetDiffuseColor(1.0000, 0.3882, 0.2784)
      self.centerlineActor.GetProperty().SetSpecularColor(1, 1, 1)
      self.centerlineActor.GetProperty().SetSpecular(.4)
      self.centerlineActor.GetProperty().SetSpecularPower(50)
      
      self.ren.AddActor(self.centerlineActor)
      self.vtkRedraw()
      
   def vtkPlotContour(self, vtkContour):
      "Affichage d'un contour"
      
      self.ren.RemoveActor(self.contourActor)
      self.vtkRedraw()
      self.contourActor = None
      
      self.contourMapper = vtk.vtkPolyDataMapper()
      if vtk.VTK_MAJOR_VERSION <= 5:
          self.contourMapper.SetInput(vtkContour.GetOutput())
      else:
          self.contourMapper.SetInputConnection(vtkContour.GetOutputPort())
      
      self.contourActor = vtk.vtkActor()
      self.contourActor.SetMapper(self.contourMapper)
      self.contourActor.GetProperty().SetDiffuseColor(0.0000, 1.0, 0.0)
      self.contourActor.GetProperty().SetSpecularColor(1, 1, 1)
      self.contourActor.GetProperty().SetSpecular(.4)
      self.contourActor.GetProperty().SetSpecularPower(50)
      
      self.ren.AddActor(self.contourActor)
      self.vtkRedraw()
   
   def vtkDeleteCenterLine(self):
      "Suppression visuelle des geometrie des centerLines"
      
      self.ren.RemoveActor(self.contourActor)
      self.vtkRedraw()
      self.contourActor = None
      
      self.ren.RemoveActor(self.centerlineActor)
      self.vtkRedraw()
      self.centerlineActor = None
