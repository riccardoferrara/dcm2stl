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
import wx, vtk
from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor
try:
   from agw import pybusyinfo as PBI
except ImportError: # if it's not there locally, try the wxPython lib.
   import wx.lib.agw.pybusyinfo as PBI

import numpy as np
import os, time
import dicomUtility
import gmshUtility
import configUtility
import vtkPipeline
import centerlineMesh
#----------------------------------------------------------------------------------------------  
class dcm2stl(wx.Frame):
   def __init__(self, id):
      #-----------------------------------------
      # Configuration de l'application
      #-----------------------------------------
      self.configUtil = configUtility.config(os.path.join(os.getcwd(),'config'))
      self.configUtil.importConfig()
      self.config = self.configUtil.config
      self.log    = configUtility.logData(os.path.join(os.getcwd(),'config'))
      self.log.logMessInfo("Demarrage de l'application.")
      #-----------------------------------------
      # Demarrage de l'application
      #-----------------------------------------
      wx.Frame.__init__(self, id, title='dcm2stl for IRRAs Technology', style =wx.DEFAULT_FRAME_STYLE, size=wx.Size(int(self.config['APP_SIZE_X']), int(self.config['APP_SIZE_Y'])))
      self.Centre()
      #-----------------------------------------
      # Add a menubar
      #-----------------------------------------
      menuBar = wx.MenuBar()
      self.openMenu = wx.Menu()
      self.viewMenu = wx.Menu()
      self.toolsMenu = wx.Menu()
      self.sliderMenu = wx.Menu()
      self.meshMenu = wx.Menu()
      #-----------------------------------------
      # Menu File
      #-----------------------------------------
      menuBar.Append(self.openMenu, "File")
      self.NewFiles  = self.openMenu.Append(wx.ID_NEW, 'N&ew\tCtrl-Shift-N', 'New reconstruction')
      self.ImportSTL = self.openMenu.Append(wx.ID_FILE, 'I&mport STL\tCtrl-Shift-S', 'Open stl files')
      self.OpenFiles = self.openMenu.Append(wx.ID_OPEN, 'I&mport dicom\tCtrl-Shift-O', 'Open dicom files')
      self.openMenu.AppendSeparator()
      self.SaveSurface     = self.openMenu.Append(wx.ID_SAVE, 'E&xport Surface\tCtrl-Shift-E', 'Export Surface with stl format')
      self.SaveSurfaceMesh = self.openMenu.Append(wx.ID_ANY, 'E&xport Mesh Surface\tCtrl-Shift-M', 'Export Mesh with stl format')
      self.openMenu.AppendSeparator()
      self.Exit    = self.openMenu.Append(wx.ID_EXIT, "E&xit\tCtrl-Shift-X", "Exit this application")
      self.Propos  = self.openMenu.Append(wx.ID_ANY, 'A& propos\tCtrl-Shift-A', 'Display message from IRRAs Technology') 
      #-----------------------------------------
      # Menu View
      #-----------------------------------------
      menuBar.Append(self.viewMenu, "View")
      self.viewAll = self.viewMenu.Append(wx.ID_ANY, 'V&iew all\tCtrl-A', 'View all objects')
      self.viewMenu.AppendSeparator()
      self.viewX    = self.viewMenu.Append(wx.ID_ANY, "V&iew X\tCtrl-X", "View in X direction")
      self.viewY    = self.viewMenu.Append(wx.ID_ANY, "V&iew Y\tCtrl-Y", "View in Y direction")
      self.viewZ    = self.viewMenu.Append(wx.ID_ANY, "V&iew Z\tCtrl-Z", "View in Z direction")
      #-----------------------------------------
      # Menu Outils
      #-----------------------------------------
      menuBar.Append(self.toolsMenu, "Data view")
      self.toolsAxis  = self.toolsMenu.Append(wx.ID_ANY, 'A&xis\tA', 'View axis', kind=wx.ITEM_CHECK)
      self.toolsBox  = self.toolsMenu.Append(wx.ID_ANY, 'B&ounding box\tB', 'View bounding box', kind=wx.ITEM_CHECK)
      self.toolsVOI  = self.toolsMenu.Append(wx.ID_ANY, 'V&OI box\tV', 'View Volum of Interest', kind=wx.ITEM_CHECK)
      self.toolsMenu.AppendSeparator()
      self.toolsSagittal  = self.toolsMenu.Append(wx.ID_ANY, 'S&agittal slice\tAlt-S', 'View sagittal slice', kind=wx.ITEM_CHECK)
      self.toolsCoronal  = self.toolsMenu.Append(wx.ID_ANY, 'C&oronal slice\tAlt-C', 'View coronal slice', kind=wx.ITEM_CHECK)
      self.toolsTransve  = self.toolsMenu.Append(wx.ID_ANY, 'T&ransversal slice\tAlt-T', 'View transversal slice', kind=wx.ITEM_CHECK)
      self.toolsMenu.AppendSeparator()
      self.toolsArtery  = self.toolsMenu.Append(wx.ID_ANY, 'A&rtery contour\tC', 'View artery contour', kind=wx.ITEM_CHECK)
      self.toolsInfoGeo = self.toolsMenu.Append(wx.ID_ANY, 'G&eometric information\tG', 'View geometric information', kind=wx.ITEM_CHECK)
      self.toolsInfoPat = self.toolsMenu.Append(wx.ID_ANY, 'P&atient information\tP', 'View patient information', kind=wx.ITEM_CHECK)
      self.toolsInfoImportedStl = self.toolsMenu.Append(wx.ID_ANY, 'I&mported .stl\tI', 'View imported .stl geometry', kind=wx.ITEM_CHECK)
      #-----------------------------------------
      # Menu Slider
      #-----------------------------------------
      menuBar.Append(self.sliderMenu, "Dicom tools")
      self.sliderSlice   = self.sliderMenu.Append(wx.ID_ANY, 'Moving medical slice', 'Sliding slice')
      self.sliderMenu.AppendSeparator()
      self.sliderVOI  = self.sliderMenu.Append(wx.ID_ANY, 'Selecting VOI', 'Selected VOI')
      self.sliderMenu.AppendSeparator()
      self.sliderArtery    = self.sliderMenu.Append(wx.ID_ANY, 'Setting isocontour', 'Sliding artery')
      self.buttonAllSurf   = self.sliderMenu.Append(wx.ID_ANY, 'G&et all surface\tH', 'Return all surface')
      self.buttonLargeSurf = self.sliderMenu.Append(wx.ID_ANY, 'G&et largest surface\tL', 'Return the largest surface')
      #-----------------------------------------
      # Menu MeshGeneration
      #-----------------------------------------
      menuBar.Append(self.meshMenu, "Mesh analysis")
      self.buttonRemeshing      = self.meshMenu.Append(wx.ID_ANY, 'Remeshing with gmsh', 'Remesh the polydata with gmsh')
      self.buttonCenterline     = self.meshMenu.Append(wx.ID_ANY, 'Compute CenterLines', 'Compute all centerline from mesh data')
      #~ self.buttonExtrudeMesh    = self.meshMenu.Append(wx.ID_ANY, 'Extrude Surface with gmsh', 'Extrude the gmsh surface')
      self.meshMenu.AppendSeparator()
      self.toolsMesh  = self.meshMenu.Append(wx.ID_ANY, 'M&esh contour\tShift-M', 'View gmsh mesh contour', kind=wx.ITEM_CHECK)
      #~ self.toolsCenterline  = self.meshMenu.Append(wx.ID_ANY, 'C&enterLine\tShift-C', 'View centerline', kind=wx.ITEM_CHECK)
      #~ self.toolsBoundary  = self.meshMenu.Append(wx.ID_ANY, 'B&oundary\tShift-B', 'View boundary', kind=wx.ITEM_CHECK)
      #~ self.toolsSolidMesh  = self.meshMenu.Append(wx.ID_ANY, 'S&olid mesh contour\tShift-S', 'View gmsh solid mesh', kind=wx.ITEM_CHECK)
      #-----------------------------------------
      # Barre de status et barre de menu pour information utilisateur
      #-----------------------------------------
      self.SetMenuBar(menuBar)
      self.sb = self.CreateStatusBar()
      self.sb.SetStatusText('Ready to work with dicom')
      #-----------------------------------------
      # Add the vtk window widget
      #-----------------------------------------
      style =  vtk.vtkInteractorStyleTrackballCamera()
      self.widget = wxVTKRenderWindowInteractor(self, -1, size = self.GetSize())
      self.widget.SetInteractorStyle(style)
      self.widget.Enable(1)
      #-----------------------------------------
      # Add a renderer
      #-----------------------------------------
      self.ren = vtk.vtkRenderer()
      self.widget.GetRenderWindow().AddRenderer(self.ren)
      self.vtkPipe = vtkPipeline.vtkPipeline(widget =  self.widget, ren = self.ren, config = self.config)
      #-----------------------------------------
      # Layout
      #-----------------------------------------
      sizer = wx.BoxSizer(wx.VERTICAL)
      sizer.Add(self.widget, 1, wx.EXPAND)
      self.SetSizer(sizer)
      self.Layout()
      #-----------------------------------------
      # Bind the menu file
      #-----------------------------------------
      self.Bind(wx.EVT_MENU, self.onNewFiles, self.NewFiles)
      self.Bind(wx.EVT_MENU, self.onOpenFiles, self.OpenFiles)
      self.Bind(wx.EVT_MENU, self.onImportSTL, self.ImportSTL)
      self.Bind(wx.EVT_MENU, self.onSaveSurface, self.SaveSurface)
      self.Bind(wx.EVT_MENU, self.onSaveSurfaceMesh, self.SaveSurfaceMesh)
      self.Bind(wx.EVT_MENU, self.onPropos, self.Propos)
      self.Bind(wx.EVT_MENU, self.onExit, self.Exit)
      #-----------------------------------------
      # Bind the menu tools
      #-----------------------------------------
      self.Bind(wx.EVT_MENU, self.onToolsSagittal, self.toolsSagittal)
      self.Bind(wx.EVT_MENU, self.onToolsCoronal, self.toolsCoronal)
      self.Bind(wx.EVT_MENU, self.onToolsTransve, self.toolsTransve)
      self.Bind(wx.EVT_MENU, self.onToolsAxis, self.toolsAxis)
      self.Bind(wx.EVT_MENU, self.onToolsBox, self.toolsBox)
      self.Bind(wx.EVT_MENU, self.onToolsVOI, self.toolsVOI)
      self.Bind(wx.EVT_MENU, self.onToolsArtery, self.toolsArtery)
      self.Bind(wx.EVT_MENU, self.onToolsInfoGeo, self.toolsInfoGeo)
      self.Bind(wx.EVT_MENU, self.onToolsInfoPat, self.toolsInfoPat)
      self.Bind(wx.EVT_MENU, self.onToolsInfoImportedStl, self.toolsInfoImportedStl)
      #-----------------------------------------
      # Bind the menu view
      #-----------------------------------------
      self.Bind(wx.EVT_MENU, self.onViewAll, self.viewAll)
      self.Bind(wx.EVT_MENU, self.onViewX, self.viewX)
      self.Bind(wx.EVT_MENU, self.onViewY, self.viewY)
      self.Bind(wx.EVT_MENU, self.onViewZ, self.viewZ)
      #-----------------------------------------
      # Bind the menu slider
      #-----------------------------------------
      self.Bind(wx.EVT_MENU, self.onSliderSlice, self.sliderSlice)
      self.Bind(wx.EVT_MENU, self.onSliderArtery, self.sliderArtery)
      self.Bind(wx.EVT_MENU, self.onSliderVOI, self.sliderVOI)
      self.Bind(wx.EVT_MENU, self.onButtonAllSurf, self.buttonAllSurf)
      self.Bind(wx.EVT_MENU, self.onButtonLargeSurf, self.buttonLargeSurf)
      #-----------------------------------------
      # Bind the menu meshing
      #-----------------------------------------
      self.Bind(wx.EVT_MENU, self.onButtonRemeshing, self.buttonRemeshing)
      self.Bind(wx.EVT_MENU, self.onButtonCenterline, self.buttonCenterline)
      #~ self.Bind(wx.EVT_MENU, self.onButtonExtrudeMesh, self.buttonExtrudeMesh)
      #~ self.Bind(wx.EVT_MENU, self.onToolsMesh, self.toolsMesh)
      #~ self.Bind(wx.EVT_MENU, self.onToolsCenterline, self.toolsCenterline)
      #~ self.Bind(wx.EVT_MENU, self.onToolsBoundary, self.toolsBoundary)
      #~ self.Bind(wx.EVT_MENU, self.onToolsSolidMesh, self.toolsSolidMesh)
      #-----------------------------------------
      # Configuration de l'application
      #-----------------------------------------
      self.initialize()
      self.menuActivate(False, 'All')
      self.dcmUtil        = dicomUtility.dicomUtility(self.config['TMP_PATH'], self.log)
      self.gmshUtil       = gmshUtility.gmshUtility(self.config['GMSH_PATH'], self.log)
      self.centerline     = centerlineMesh.centerlineMesh(self.config['GMSH_PATH'], self.log)
   
   def menuActivate(self, boolean, type):
      "Active ou desactive les menu "
      #-----------------------------------------
      # Menu Outils
      #-----------------------------------------
      if type == 'Outils' or type == 'All' :
         self.toolsAxis.Enable(boolean)
         self.toolsBox.Enable(boolean)
         self.toolsVOI.Enable(boolean)
         self.toolsSagittal.Enable(boolean)
         self.toolsCoronal.Enable(boolean)
         self.toolsTransve.Enable(boolean)
         self.toolsArtery.Enable(boolean)
         self.toolsInfoGeo.Enable(boolean)
         self.toolsInfoPat.Enable(boolean)
      if type == 'ImportedStl' or type == 'All' :
         self.toolsInfoImportedStl.Enable(boolean)
      #-----------------------------------------
      # Menu Slider
      #-----------------------------------------
      if type == 'Slider' or type == 'All' :
         self.sliderSlice.Enable(boolean)
         self.sliderVOI.Enable(boolean)
         self.sliderArtery.Enable(boolean)
         self.buttonAllSurf.Enable(boolean)
         self.buttonLargeSurf.Enable(boolean)
      #-----------------------------------------
      # Menu MeshGeneration
      #-----------------------------------------
      if type == 'Mesh' or type == 'All' :
         self.buttonRemeshing.Enable(boolean)
         self.buttonCenterline.Enable(boolean)
         self.toolsMesh.Enable(boolean)
         #~ self.toolsSolidMesh.Enable(boolean)
      #-----------------------------------------
      # Menu View
      #-----------------------------------------
      if type == 'View' or type == 'All' :
         self.viewAll.Enable(boolean)
         self.viewX.Enable(boolean)
         self.viewY.Enable(boolean)
         self.viewZ.Enable(boolean)
         
   def onButtonRemeshing(self, evt):
      " Remeshing with gmsh application "
      if self.toolsArtery.IsChecked():
         # Si le contour existe
         self.PleaseCompute(evt)
         self.vtkPipe.vtkWriteSTLfile(os.path.join(self.config['GMSH_PATH'],'Surface.stl'), self.config['STL_TYPE'],'Surface')
         self.gmshUtil.Remeshing()
         self.vtkPipe.vtkReadVTKfile(os.path.join(self.config['GMSH_PATH'],'SurfaceMesh.vtk'))
         self.vtkPipe.vtkPlotMesh()
         #~ self.vtkPipe.vtkImportSTLfile(os.path.join(self.config['GMSH_PATH'],'SurfaceMesh.stl'))
         self.meshMenu.Check(self.toolsMesh.GetId(), True)
         self.DestroyCompute()
      else : 
         # Si le contour n'existe pas
         text = "WARNING : Please select the surface to extrude"
         self.log.logMessWarm(text)
         self.onProblem(text)
   
   def onButtonExtrudeMesh(self, evt):
      " Extrude with gmsh application "
      if self.toolsMesh.IsChecked():
         # Si le contour existe
         self.PleaseCompute(evt)
         #~ self.vtkPipe.vtkWriteSTLfile(os.path.join(self.config['GMSH_PATH'],'Surface.stl'), self.config['STL_TYPE'],'Surface')
         self.gmshUtil.Extrude()
         self.vtkPipe.vtkReadVTKSolidfile(os.path.join(self.config['GMSH_PATH'],'SolidMesh.vtk'))
         self.meshMenu.Check(self.toolsSolidMesh.GetId(), True)
         self.DestroyCompute()
      else : 
         # Si le contour n'existe pas
         text = "WARNING : Please select the surface to extrude"
         self.log.logMessWarm(text)
         self.onProblem(text)

   def onExit(self, evt):
      """Fermeture du programme."""
      self.dcmUtil.remove_Images()
      self.gmshUtil.remove_Mesh()
      
      self.config['VALSEUIL'] = str(self.vtkPipe.valSeuil)
      self.config['VALFILTER'] = str(self.vtkPipe.valFilter)
      self.config['VALGAUSSIANFILTER'] = str(self.vtkPipe.valGaussianFilter)
      self.config['VALSMOOTHITER'] = str(self.vtkPipe.valSmoothIter)
      self.config['VALDECIMATE'] = str(int(self.vtkPipe.valDecimate*100))
      self.config['APP_SIZE_X'] = str(self.GetSize()[0])
      self.config['APP_SIZE_Y'] = str(self.GetSize()[1])
      
      self.configUtil.saveConfig(self.config)
      self.Close()

   def initialize(self):
      """ Initialisation de l'application """
      
      self.Description = """dicom2stl is an advanced dicom viewer for 
      the Unix operating system. Features include powerful built-in surface 
      reconstruction from medical image IRM or CT.
      """
      self.Licence = """This software is governed by the CeCILL license under French law and 
      abiding by the rules of distribution of free software.  You can  use, 
      modify and/ or redistribute the software under the terms of the CeCILL 
      license as circulated by CEA, CNRS and INRIA at the following URL 
      "http://www.cecill.info"."""
      
      self.info = wx.AboutDialogInfo()
      self.info.SetIcon(wx.Icon(os.path.join(os.getcwd(),os.path.join('ico',os.path.join('128x128','ViewResults.png'))), wx.BITMAP_TYPE_PNG))
      self.info.SetName('dicom2stl')
      self.info.SetVersion(self.configUtil.version)
      self.info.SetDescription(self.Description)
      self.info.SetCopyright('(C) 2014 - 2015 dominique ambard')
      self.info.SetWebSite('http://www.lmgc.univ-montp2.fr/')
      self.info.SetLicence(self.Licence)
      self.info.AddDeveloper('dominique ambard')
      self.info.AddDocWriter('dominique ambard')
      self.info.AddArtist('dominique ambard')
      self.info.AddTranslator('dominique ambard')
      self.Compute = None

   def onPropos(self, evt):
      """ Affiche une fenetre d'information sur le logiciel"""
      wx.AboutBox(self.info)

   def onProblem(self, msg):
      """ Affichage d'un probleme """
      dlg = wx.MessageDialog(self, msg,
                               'A problem is detect',
                               wx.OK | wx.ICON_INFORMATION
                               #wx.YES_NO | wx.NO_DEFAULT | wx.CANCEL | wx.ICON_INFORMATION
                               )
      dlg.ShowModal()
      dlg.Destroy()

   def startProgressBar(self, func, *args):
      """ Demarre une barre d'avancement du chargement des images"""
      thread = threading.Thread(target = func, args = args)
      thread.setDaemon(True)
      thread.start()

   def onNewFiles(self, evt):
      """ Nouvelle reconstruction """
      self.dcmUtil.remove_Images()
      self.gmshUtil.remove_Mesh()
      self.vtkPipe.vtkDeleteIsoContour()
      self.vtkPipe.vtkDeleteBlockFilterIsoContour()
      self.vtkPipe.vtkDeletePlotIsoContour()
      self.vtkPipe.vtkDeleteAllActor()
      self.vtkPipe.vtkRedraw()
      self.menuActivate(False, 'All')

   def onOpenFiles(self, evt):
      """ Chargement des fichiers dicom"""
      wildcard = "Dicom source (*.dcm)|*.dcm|"     \
                 "Dicom source (*.DCM)|*.DCM|" \
                 "All files (*)|*"
      self.path = []
      self.listFiles = []
      self.dcmUtil.remove_Images()
      self.gmshUtil.remove_Mesh()
      #-----------------------------------------
      # Interface de selection des fichiers
      #-----------------------------------------
      dlg = wx.FileDialog(self, "Choose a list of file", os.getcwd(), " "*80, wildcard, wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR)
      
      if dlg.ShowModal() == wx.ID_OK:
         path = dlg.GetPath()
         self.path = os.path.dirname(path)
         self.listFiles = dlg.GetFilenames()
         self.sb.SetStatusText("You selected: %s" % self.path)
      dlg.Destroy()
      #-----------------------------------------
      # Progression du chargement
      #-----------------------------------------
      self.import_all_dicom()
      self.menuActivate(True, 'All')
      self.vtkPipe.vtkInitialData(self.config['TMP_PATH'], self.config['SURFACE_TYPE'])
      self.vtkPipe.vtkPlotBoundingBox()
      self.toolsMenu.Check(self.toolsBox.GetId(), True)
      self.vtkPipe.vtkInfoGeometric()
      self.toolsMenu.Check(self.toolsInfoGeo.GetId(), True)
      
   def onToolsSagittal(self, evt):
      """ Affichage de la coupe sagittale"""
      #-----------------------------------------
      # Creation d'une coupe sagittale
      #-----------------------------------------
      if self.toolsSagittal.IsChecked():
         self.vtkPipe.vtkSagittalSlice()
      else:
         self.vtkPipe.vtkDeleteSagittalSlice()

   def onToolsCoronal(self, evt):
      """ Affichage de la coupe coronale"""
      #-----------------------------------------
      # Creation d'une coupe coronale
      #-----------------------------------------
      if self.toolsCoronal.IsChecked():
         self.vtkPipe.vtkCoronalSlice()
      else:
         self.vtkPipe.vtkDeleteCoronalSlice()

   def onToolsTransve(self, evt):
      """ Affichage de la coupe transversale"""
      #-----------------------------------------
      # Creation d'une coupe transversale
      #-----------------------------------------
      if self.toolsTransve.IsChecked():
         self.vtkPipe.vtkTransversalSlice()
      else:
         self.vtkPipe.vtkDeleteTransversalSlice()

   def onToolsAxis(self, evt):
      """ Affichage du repere"""
      #-----------------------------------------
      # Creation d'un repere de travail
      #-----------------------------------------
      if self.toolsAxis.IsChecked():
         self.vtkPipe.vtkAxis()
      else:
         self.vtkPipe.vtkDeleteAxis()
         
   def onToolsVOI(self, evt):
      """ Affichage de la zone de travail"""
      #-----------------------------------------
      # Creation de la bounding box de la VOI
      #-----------------------------------------
      if self.toolsVOI.IsChecked():
         self.vtkPipe.vtkPlotBoundingVOI()
      else:
         self.vtkPipe.vtkDeleteBoundingVOI()
   
   def onToolsBox(self, evt):
      """ Affichage de la zone de travail"""
      #-----------------------------------------
      # Creation de la bounding box du volume total
      #-----------------------------------------
      if self.toolsBox.IsChecked():
         self.vtkPipe.vtkPlotBoundingBox()
      else:
         self.vtkPipe.vtkDeleteBoundingBox()
         
   def onToolsMesh(self, evt):
      """ Affichage du maillage gmsh """
      #-----------------------------------------
      # Creation du maillage EF des arteres
      #-----------------------------------------
      if self.toolsMesh.IsChecked():
         self.vtkPipe.vtkPlotMesh()
      else:
         self.vtkPipe.vtkDeleteMesh()

   def onToolsSolidMesh(self, evt):
      """ Affichage du maillage 3D gmsh """
      #-----------------------------------------
      # Creation du maillage EF des arteres
      #-----------------------------------------
      if self.toolsSolidMesh.IsChecked():
         self.vtkPipe.vtkReadVTKSolidfile(os.path.join(self.config['GMSH_PATH'],'SolidMesh.vtk'))
      else:
         self.vtkPipe.vtkDeleteSolidMesh()

   def onToolsArtery(self, evt):
      """ Affichage de la zone de travail"""
      #-----------------------------------------
      # Creation du contour des arteres
      #-----------------------------------------
      if self.toolsArtery.IsChecked():
         self.PleaseCompute(evt)
         wx.SetCursor(wx.StockCursor(wx.CURSOR_WAIT))
         #-----------------------------------------
         # Creation de la donnee surface
         #-----------------------------------------
         self.vtkPipe.vtkCreateIsoContour(config = self.config['SEGMENTATION_TYPE'])
         #-----------------------------------------
         # Creation du filtrage de la surface
         #-----------------------------------------
         self.vtkPipe.vtkCreateFilterIsoContour()
         #-----------------------------------------
         # Creation des infos geometrique
         #-----------------------------------------
         self.vtkPipe.vtkComputeGeometricSurface()
         self.vtkPipe.vtkCreateGeometricInformations()
         #-----------------------------------------
         # Affichage du contour surfacique
         #-----------------------------------------
         self.vtkPipe.vtkPlotIsoContour()
         self.vtkPipe.vtkRedraw()
         wx.SetCursor(wx.StockCursor(wx.CURSOR_ARROW))
         self.DestroyCompute()
      else:
         #-----------------------------------------
         # Effacement du contour
         #-----------------------------------------
         self.vtkPipe.vtkDeleteIsoContour()
         self.vtkPipe.vtkDeleteBlockFilterIsoContour()
         self.vtkPipe.vtkDeletePlotIsoContour()

   def onToolsInfoGeo(self, evt):
      "Affichage des informations geometriques"
      if self.toolsInfoGeo.IsChecked():
         self.vtkPipe.vtkInfoGeometric()
      else :
         self.vtkPipe.vtkDeleteInfoGeometric()

   def onToolsInfoPat(self, evt):
      "Affichage des informations patients"
      if self.toolsInfoPat.IsChecked():
         self.vtkPipe.vtkInfoPatient()
      else :
         self.vtkPipe.vtkDeleteInfoPatient()
         
   def onToolsInfoImportedStl(self, evt):
      "Affichage de la geometrie .stl importée"
      if self.toolsInfoImportedStl.IsChecked():
         self.vtkPipe.vtkImportSTLfile(self.stlPath)
#         self.vtkPipe.valOpacity = 2
#         self.vtkPipe.meshActor.GetProperty().SetOpacity(self.vtkPipe.valOpacity/10.0)
#         self.vtkPipe.vtkRedraw()
      else :
         self.vtkPipe.vtkDeleteImportedStl()
         
         
   def onViewAll(self, evt):
      "Gestion de l'angle de visualisation"
      self.vtkPipe.ren.ResetCamera()
      self.vtkPipe.vtkRedraw()

   def onViewZ(self, evt):
      "Gestion de l'angle de visualisation"
      self.vtkPipe.vtkCameraViewZ()

   def onViewY(self, evt):
      "Gestion de l'angle de visualisation"
      self.vtkPipe.vtkCameraViewY()

   def onViewX(self, evt):
      "Gestion de l'angle de visualisation"
      self.vtkPipe.vtkCameraViewX()
      
   def import_all_dicom(self):
      "Lecture du fichier dicom"
      #-----------------------------------------
      # Boucle sur les fichiers dicom
      #-----------------------------------------
      self.vtkPipe.nb_coupe = len(self.listFiles)
      
      ProgressLoading = wx.ProgressDialog('IRRAs Technology  Import dicom files', 'Please wait...' + ' '*40 ,maximum = self.vtkPipe.nb_coupe, parent=self,
                               style = wx.PD_APP_MODAL
                                     | wx.PD_ELAPSED_TIME
                                    #| wx.PD_ESTIMATED_TIME
                                    #| wx.PD_REMAINING_TIME
                                          )
      n = 0
      keepGoing = True
      
      while keepGoing and n < self.vtkPipe.nb_coupe:

         data = self.dcmUtil.import_dicom(self.path, self.listFiles[n])
         #-----------------------------------------
         if n==0 : 
            # Allocation de la visu vtk
            self.vtkPipe.lig,self.vtkPipe.col = self.dcmUtil.get_slice_size(data)
            self.vtkPipe.pixel_spacing_x, self.vtkPipe.pixel_spacing_y = self.dcmUtil.get_pixel_spacing(data)
            self.vtkPipe.slice_spacing = self.dcmUtil.get_slice_spacing(data)
            position      = self.dcmUtil.get_slice_position(data)
            indice_coupe  = self.dcmUtil.get_slice_number(data)
            self.vtkPipe.min_pixel_value = self.dcmUtil.get_min_pixel_value(data)
            self.vtkPipe.max_pixel_value = self.dcmUtil.get_max_pixel_value(data)
            self.dcmUtil.add_slice2vtk(self.dcmUtil.get_slice_data(data), indice_coupe)
            self.dcmUtil.createPatientInformations(data)
            self.vtkPipe.PatientInfos = self.dcmUtil.PatientInfos
         #-----------------------------------------
         if n > 0 :
            # Recuperation d'une coupe
            indice_coupe  = self.dcmUtil.get_slice_number(data)
            self.dcmUtil.add_slice2vtk(self.dcmUtil.get_slice_data(data), indice_coupe)
            min_pixel_value = self.dcmUtil.get_min_pixel_value(data)
            max_pixel_value = self.dcmUtil.get_max_pixel_value(data)
            if self.vtkPipe.min_pixel_value > min_pixel_value : self.vtkPipe.min_pixel_value = min_pixel_value
            if self.vtkPipe.max_pixel_value < max_pixel_value : self.vtkPipe.max_pixel_value = max_pixel_value
         #-----------------------------------------
         (keepGoing, skip) = ProgressLoading.Update(n)
         print "Recuperation image  : ",n," - ",self.listFiles[n]
         n += 1
         time.sleep(0.001)
      #-----------------------------------------
      ProgressLoading.Destroy()
      
   def onSaveSurface(self, evt):
      """ Sauvegarde au format stl de la surface """
      
      if self.vtkPipe.smoothAneurysm <> None :
         wildcard = "Stl file (*.stl)|*.stl|"     \
                    "All files (*.*)|*.*"
   
         dlg = wx.FileDialog(self, message="Save file as ...", defaultDir=os.getcwd(), 
                             defaultFile="", wildcard=wildcard, style=wx.SAVE )
   
         if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            filename = dlg.GetFilename()
            #-----------------------------------------
            # Ecriture de l isocontour
            #-----------------------------------------
            self.vtkPipe.vtkWriteSTLfile(path, self.config['STL_TYPE'], 'Surface' )
         dlg.Destroy()
      else :
         text = "WARNING : Please create a surface."
         self.onProblem(text)
         self.log.logMessWarm(text)
      
   def onSaveSurfaceMesh(self, evt):
      """ Sauvegarde au format stl du maillage """
      
      if self.vtkPipe.outputMesh <> None :
         wildcard = "Stl file (*.stl)|*.stl|"     \
                    "All files (*.*)|*.*"
   
         dlg = wx.FileDialog(self, message="Save file as ...", defaultDir=os.getcwd(), 
                             defaultFile="", wildcard=wildcard, style=wx.SAVE )
   
         if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            filename = dlg.GetFilename()
            #-----------------------------------------
            # Ecriture de l isocontour
            #-----------------------------------------
            self.vtkPipe.vtkWriteSTLfile(path, self.config['STL_TYPE'],'Mesh' )
         dlg.Destroy()
      else :
         text = "WARNING : Please create a mesh."
         self.onProblem(text)
         self.log.logMessWarm(text)
      
   def onImportSTL(self, evt):
      """ Sauvegarde au format stl """
      
      wildcard = "Stl file (*.stl)|*.stl|"     \
                 "All files (*.*)|*.*"

      dlg = wx.FileDialog(self, message="Import file as ...", defaultDir=os.getcwd(), 
                        defaultFile="", wildcard=wildcard, style=wx.FD_OPEN )

      if dlg.ShowModal() == wx.ID_OK:
         stlPath = dlg.GetPath()
         self.filename = dlg.GetFilename()
         self.stlPath = os.path.dirname(stlPath)+"/"+self.filename
         self.sb.SetStatusText("You selected: %s" % self.stlPath)    
         #-----------------------------------------
         # Importation d'un fichier stl
         #-----------------------------------------
         self.vtkPipe.vtkImportSTLfile(stlPath)
         #activate the menu and check the option
         self.menuActivate(True, 'ImportedStl')
         self.toolsMenu.Check(self.toolsInfoImportedStl.GetId(), True)
      dlg.Destroy()
      
   def onSliderSlice(self, evt):
      """Affiche une boite de dialogue pour reglage par slider des coupes """
      self.sliderSliceDialog = wx.Dialog(self, id=wx.ID_ANY, title="Slice parameters ", size=wx.Size(240, 360), style=wx.DEFAULT_DIALOG_STYLE)
      #-----------------------------------------
      # Slider coronal
      #-----------------------------------------
      text = wx.StaticText(self.sliderSliceDialog, id=wx.ID_ANY, label="Sagittal slice", pos=(20, 20))
      slider0 = wx.Slider(self.sliderSliceDialog, id=wx.ID_ANY, value=self.vtkPipe.valSagittal, minValue=0, maxValue=self.vtkPipe.col-1, pos=(20, 40), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider0.Bind(wx.EVT_SCROLL, self.onSliderSagittal)
      wx.StaticLine(self.sliderSliceDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 80))
      #-----------------------------------------
      # Slider sagittal
      #-----------------------------------------
      text = wx.StaticText(self.sliderSliceDialog, id=wx.ID_ANY, label="Coronal slice", pos=(20, 100))
      slider1 = wx.Slider(self.sliderSliceDialog, id=wx.ID_ANY, value=self.vtkPipe.valCoronal, minValue=0, maxValue=self.vtkPipe.lig-1, pos=(20, 120), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider1.Bind(wx.EVT_SCROLL, self.onSliderCoronal)
      wx.StaticLine(self.sliderSliceDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 160))
      #-----------------------------------------
      # Slider transversal
      #-----------------------------------------
      text = wx.StaticText(self.sliderSliceDialog, id=wx.ID_ANY, label="Transversal slice", pos=(20, 180))
      slider2 = wx.Slider(self.sliderSliceDialog, id=wx.ID_ANY, value=self.vtkPipe.valTransversal, minValue=0, maxValue=self.vtkPipe.nb_coupe-1, pos=(20, 200), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider2.Bind(wx.EVT_SCROLL, self.onSliderTransversal)
      wx.StaticLine(self.sliderSliceDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20,240))
      
      #-----------------------------------------
      # Bouton de fermeture de la boite de dialogue
      #-----------------------------------------
      close = wx.Button(self.sliderSliceDialog, label='Close', pos=(130, 280), size=(70, 25))
      close.Bind(wx.EVT_BUTTON, self.onCloseSliderSlice)
      self.sliderSliceDialog.Show()
      
   def onCloseSliderSlice(self, evt):
      """ Fermeture de la boite de dialogue SliderArtery"""
      self.sliderSliceDialog.Close(True)  
         
   def onSliderCoronal(self, evt):
      """Modifie la coupe coronale"""
      if self.vtkPipe.coronalCut<>None:
         obj = evt.GetEventObject()
         self.vtkPipe.valCoronal = obj.GetValue()
         self.vtkPipe.coronalCut.SetDisplayExtent(int(self.vtkPipe.valCoronal), int(self.vtkPipe.valCoronal), 0, self.vtkPipe.lig-1, 0, self.vtkPipe.nb_coupe - 1)
         self.vtkPipe.vtkRedraw()
      
   def onSliderSagittal(self, evt):
      """Modifie la coupe coronale"""
      if self.vtkPipe.sagittalCut<>None:
         obj = evt.GetEventObject()
         self.vtkPipe.valSagittal = obj.GetValue()
         self.vtkPipe.sagittalCut.SetDisplayExtent(0, self.vtkPipe.col-1, int(self.vtkPipe.valSagittal), int(self.vtkPipe.valSagittal), 0, self.vtkPipe.nb_coupe - 1)
         self.vtkPipe.vtkRedraw()

   def onSliderTransversal(self, evt):
      """Modifie la coupe coronale"""
      if self.vtkPipe.transversalCut<>None:
         obj = evt.GetEventObject()
         self.vtkPipe.valTransversal = obj.GetValue()
         self.vtkPipe.transversalCut.SetDisplayExtent(0, self.vtkPipe.col-1, 0, self.vtkPipe.lig-1, int(self.vtkPipe.valTransversal), int(self.vtkPipe.valTransversal))
         self.vtkPipe.vtkRedraw()

   def onSliderArtery(self, evt):
      """Affiche une boite de dialogue pour reglage par slider du contour """
      
      if not(self.toolsArtery.IsChecked()):
         self.PleaseCompute(evt)
         wx.SetCursor(wx.StockCursor(wx.CURSOR_WAIT))
         #-----------------------------------------
         # Creation de la donnee surface
         #-----------------------------------------
         self.vtkPipe.vtkCreateIsoContour(config = self.config['SEGMENTATION_TYPE'])
         #-----------------------------------------
         # Creation du filtrage de la surface
         #-----------------------------------------
         self.vtkPipe.vtkCreateFilterIsoContour()
         #-----------------------------------------
         # Creation des infos geometrique
         #-----------------------------------------
         self.vtkPipe.vtkComputeGeometricSurface()
         self.vtkPipe.vtkCreateGeometricInformations()
         #-----------------------------------------
         # Affichage du contour surfacique
         #-----------------------------------------
         self.vtkPipe.vtkPlotIsoContour()
         self.vtkPipe.vtkRedraw()
         wx.SetCursor(wx.StockCursor(wx.CURSOR_ARROW))
         self.DestroyCompute()
         self.toolsMenu.Check(self.toolsArtery.GetId(), True)
         
      self.sliderArteryDialog = wx.Dialog(self, id=wx.ID_ANY, title="Surfacic contour parameters ", size=wx.Size(240, 600), style=wx.DEFAULT_DIALOG_STYLE)
      #-----------------------------------------
      # Slider valeur du seuil
      #-----------------------------------------
      text = wx.StaticText(self.sliderArteryDialog, id=wx.ID_ANY, label="Intensity value", pos=(20, 20))
      slider0 = wx.Slider(self.sliderArteryDialog, id=wx.ID_ANY, value=self.vtkPipe.valSeuil, minValue=int(self.vtkPipe.vtkVolumRange[0]) , maxValue=int(self.vtkPipe.vtkVolumRange[1]) , pos=(20, 40), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider0.Bind(wx.EVT_SCROLL, self.onSliderSeuil)
      wx.StaticLine(self.sliderArteryDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 80))
      #-----------------------------------------
      # Slider valeur du filtre sur l'image
      #-----------------------------------------
      text = wx.StaticText(self.sliderArteryDialog, id=wx.ID_ANY, label="Image filter level", pos=(20, 100))
      slider1 = wx.Slider(self.sliderArteryDialog, id=wx.ID_ANY, value=self.vtkPipe.valGaussianFilter, minValue=0, maxValue=4, pos=(20, 120), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider1.Bind(wx.EVT_SCROLL, self.onSliderImageFilter)
      wx.StaticLine(self.sliderArteryDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 160))
      #-----------------------------------------
      # Slider valeur du filtre sur la surface
      #-----------------------------------------
      text = wx.StaticText(self.sliderArteryDialog, id=wx.ID_ANY, label="Surface filter level", pos=(20, 180))
      slider2 = wx.Slider(self.sliderArteryDialog, id=wx.ID_ANY, value=self.vtkPipe.valFilter, minValue=0, maxValue=3, pos=(20, 200), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider2.Bind(wx.EVT_SCROLL, self.onSliderSurfFilter)
      wx.StaticLine(self.sliderArteryDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 240))
      #-----------------------------------------
      # Slider smooth surface
      #-----------------------------------------
      text = wx.StaticText(self.sliderArteryDialog, id=wx.ID_ANY, label="Smooth surface", pos=(20, 260))
      slider3 = wx.Slider(self.sliderArteryDialog, id=wx.ID_ANY, value=self.vtkPipe.valSmoothIter, minValue=0, maxValue=100, pos=(20, 280), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider3.Bind(wx.EVT_SCROLL, self.onSliderSmoothIter)
      wx.StaticLine(self.sliderArteryDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 320))
      #-----------------------------------------
      # Slider decimate
      #-----------------------------------------
      text = wx.StaticText(self.sliderArteryDialog, id=wx.ID_ANY, label="Decimate", pos=(20, 340))
      slider4 = wx.Slider(self.sliderArteryDialog, id=wx.ID_ANY, value=self.vtkPipe.valDecimate*100, minValue=0, maxValue=100, pos=(20, 360), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider4.Bind(wx.EVT_SCROLL, self.onSliderDecimate)
      wx.StaticLine(self.sliderArteryDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 400))
      #-----------------------------------------
      # Slider de l'opacite
      #-----------------------------------------
      text = wx.StaticText(self.sliderArteryDialog, id=wx.ID_ANY, label="Opacity level", pos=(20, 420))
      slider5 = wx.Slider(self.sliderArteryDialog, id=wx.ID_ANY, value=self.vtkPipe.valOpacity, minValue=0, maxValue=10, pos=(20, 440), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider5.Bind(wx.EVT_SCROLL, self.onSliderOpacity)
      wx.StaticLine(self.sliderArteryDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 480))
      #-----------------------------------------
      # Bouton de fermeture de la boite de dialogue
      #-----------------------------------------
      close = wx.Button(self.sliderArteryDialog, label='Close', pos=(130, 520), size=(70, 25))
      close.Bind(wx.EVT_BUTTON, self.onCloseSliderArtery)
      #-----------------------------------------
      # Bouton de calcul du contour
      #-----------------------------------------
      compute = wx.Button(self.sliderArteryDialog, label='Compute', pos=(40, 520), size=(70, 25))
      compute.Bind(wx.EVT_BUTTON, self.onComputeArtery)
      self.sliderArteryDialog.Show()

   def onCloseSliderArtery(self, evt):
      """ Fermeture de la boite de dialogue SliderArtery"""
      self.sliderArteryDialog.Close(True) 
      
   def onComputeArtery(self,evt):
      """ Applique les modifications sur le contour """
      self.PleaseCompute(evt)
      self.vtkPipe.vtkComputeGeometricSurface()
      self.vtkPipe.vtkCreateGeometricInformations()
      self.vtkPipe.vtkCreateFilterIsoContour()
      if self.toolsInfoGeo.IsChecked() : 
         self.vtkPipe.vtkInfoGeo.SetInput(self.vtkPipe.GeoInfos)
      self.vtkPipe.vtkRedraw()
      self.DestroyCompute()
      self.vtkPipe.vtkHistogram()

   def onSliderSeuil(self, evt):
      """Modifie le seuil de la creation du contour"""
      if self.vtkPipe.aneurysm<>None:
         obj = evt.GetEventObject()
         self.vtkPipe.valSeuil = obj.GetValue()
         #-----------------------------------------
         # Modification du seuil de l isocontour
         #-----------------------------------------
         self.vtkPipe.aneurysmExtractor.SetValue(0, self.vtkPipe.valSeuil)
      
   def onSliderSurfFilter(self, evt):
      """Modifie la taille de la discretisation surfacique"""
      if self.vtkPipe.aneurysm<>None:
         obj = evt.GetEventObject()
         self.vtkPipe.valFilter = obj.GetValue()
         #-----------------------------------------
         # Modification du filtre
         #-----------------------------------------
         self.vtkPipe.subdiviseAneurysm.SetNumberOfSubdivisions(self.vtkPipe.valFilter)

   def onSliderImageFilter(self, evt):
      """Modifie la valeur du flou gaussien de l'image"""
      if self.vtkPipe.aneurysm<>None:
         obj = evt.GetEventObject()
         self.vtkPipe.valGaussianFilter = obj.GetValue()
         #-----------------------------------------
         # Modification du filtre
         #-----------------------------------------
         self.vtkPipe.vtkVolumBlur.SetStandardDeviation(self.vtkPipe.valGaussianFilter)

   def onSliderDecimate(self, evt):
      """Modifie le seuil de la creation du contour"""
      if self.vtkPipe.aneurysm<>None:
         obj = evt.GetEventObject()
         self.vtkPipe.valDecimate = obj.GetValue()/100.0
         #-----------------------------------------
         # Modification de decimate
         #-----------------------------------------
         self.vtkPipe.aneurysmDecimator.SetTargetReduction(self.vtkPipe.valDecimate) 
         
   def onSliderSmoothIter(self, evt):
      """Modifie le seuil de la creation du contour"""
      if self.vtkPipe.aneurysm<>None:
         obj = evt.GetEventObject()
         self.vtkPipe.valSmoothIter = obj.GetValue()
         #-----------------------------------------
         # Modification du smooth surface
         #-----------------------------------------
         self.vtkPipe.smoothAneurysm.SetNumberOfIterations(self.vtkPipe.valSmoothIter)

   def onSliderOpacity(self, evt):
      """Gere l'opacite de la construction surfacique"""
      obj = evt.GetEventObject()
      self.vtkPipe.valOpacity = obj.GetValue()
      self.vtkPipe.aneurysm.GetProperty().SetOpacity(self.vtkPipe.valOpacity/10.0)
      self.vtkPipe.vtkRedraw()

   def onSliderVOI(self, evt):
      """Affiche une boite de dialogue pour reglage par slider du VOI """
      
      if not(self.toolsVOI.IsChecked()):
         self.vtkPipe.vtkPlotBoundingVOI()
         self.toolsMenu.Check(self.toolsVOI.GetId(), True)
      
      self.sliderVOIDialog = wx.Dialog(self, id=wx.ID_ANY, title="VOI parameters ", size=wx.Size(240, 600), style=wx.DEFAULT_DIALOG_STYLE)
      #-----------------------------------------
      # Slider de la position du point 1 sur X
      #-----------------------------------------
      text = wx.StaticText(self.sliderVOIDialog, id=wx.ID_ANY, label="Point 1 - X position", pos=(20, 20))
      slider0 = wx.Slider(self.sliderVOIDialog, id=wx.ID_ANY, value=self.vtkPipe.VOI_pt1[0], minValue=0 , maxValue=self.vtkPipe.col , pos=(20, 40), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider0.Bind(wx.EVT_SCROLL, self.onSliderVOI_Pt1_X)
      #-----------------------------------------
      # Slider de la position du point 1 sur Y
      #-----------------------------------------
      text = wx.StaticText(self.sliderVOIDialog, id=wx.ID_ANY, label="Point 1 - Y position", pos=(20, 100))
      slider1 = wx.Slider(self.sliderVOIDialog, id=wx.ID_ANY, value=self.vtkPipe.VOI_pt1[1], minValue=0, maxValue=self.vtkPipe.lig, pos=(20, 120), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider1.Bind(wx.EVT_SCROLL, self.onSliderVOI_Pt1_Y)
      #-----------------------------------------
      # Slider de la position du point 1 sur Z
      #-----------------------------------------
      text = wx.StaticText(self.sliderVOIDialog, id=wx.ID_ANY, label="Point 1 - Z position", pos=(20, 180))
      slider2 = wx.Slider(self.sliderVOIDialog, id=wx.ID_ANY, value=self.vtkPipe.VOI_pt1[2], minValue=0, maxValue=self.vtkPipe.nb_coupe, pos=(20, 200), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider2.Bind(wx.EVT_SCROLL, self.onSliderVOI_Pt1_Z)
      wx.StaticLine(self.sliderVOIDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 260))
      #-----------------------------------------
      # Slider de la position du point 2 sur X
      #-----------------------------------------
      text = wx.StaticText(self.sliderVOIDialog, id=wx.ID_ANY, label="Point 2 - X position", pos=(20, 280))
      slider3 = wx.Slider(self.sliderVOIDialog, id=wx.ID_ANY, value=self.vtkPipe.VOI_pt2[0], minValue=0 , maxValue=self.vtkPipe.col , pos=(20, 300), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider3.Bind(wx.EVT_SCROLL, self.onSliderVOI_Pt2_X)
      #-----------------------------------------
      # Slider de la position du point 2 sur Y
      #-----------------------------------------
      text = wx.StaticText(self.sliderVOIDialog, id=wx.ID_ANY, label="Point 2 - Y position", pos=(20, 360))
      slider4 = wx.Slider(self.sliderVOIDialog, id=wx.ID_ANY, value=self.vtkPipe.VOI_pt2[1], minValue=0, maxValue=self.vtkPipe.lig, pos=(20, 380), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider4.Bind(wx.EVT_SCROLL, self.onSliderVOI_Pt2_Y)
      #-----------------------------------------
      # Slider de la position du point 2 sur Z
      #-----------------------------------------
      text = wx.StaticText(self.sliderVOIDialog, id=wx.ID_ANY, label="Point 2 - Z position", pos=(20, 440))
      slider5 = wx.Slider(self.sliderVOIDialog, id=wx.ID_ANY, value=self.vtkPipe.VOI_pt2[2], minValue=0, maxValue=self.vtkPipe.nb_coupe, pos=(20, 460), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider5.Bind(wx.EVT_SCROLL, self.onSliderVOI_Pt2_Z)
      wx.StaticLine(self.sliderVOIDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 520))
      #-----------------------------------------
      # Bouton de fermeture de la boite de dialogue
      #-----------------------------------------
      close = wx.Button(self.sliderVOIDialog, label='Close', pos=(130, 540), size=(70, 25))
      close.Bind(wx.EVT_BUTTON, self.onCloseSliderVOI)
      self.sliderVOIDialog.Show()

   def onSliderVOI_Pt1_X(self, evt):
      """ Fonction pour modifier la position de la VOI"""
      obj = evt.GetEventObject()
      self.vtkPipe.VOI_pt1[0] = obj.GetValue()
      self.vtkPipe.vtkMoveVOI()
      
   def onSliderVOI_Pt1_Y(self, evt):
      """ Fonction pour modifier la position de la VOI"""
      obj = evt.GetEventObject()
      self.vtkPipe.VOI_pt1[1] = obj.GetValue()
      self.vtkPipe.vtkMoveVOI()
      
   def onSliderVOI_Pt1_Z(self, evt):
      """ Fonction pour modifier la position de la VOI"""
      obj = evt.GetEventObject()
      self.vtkPipe.VOI_pt1[2] = obj.GetValue()
      self.vtkPipe.vtkMoveVOI()
      
   def onSliderVOI_Pt2_X(self, evt):
      """ Fonction pour modifier la position de la VOI"""
      obj = evt.GetEventObject()
      self.vtkPipe.VOI_pt2[0] = obj.GetValue()
      self.vtkPipe.vtkMoveVOI()
      
   def onSliderVOI_Pt2_Y(self, evt):
      """ Fonction pour modifier la position de la VOI"""
      obj = evt.GetEventObject()
      self.vtkPipe.VOI_pt2[1] = obj.GetValue()
      self.vtkPipe.vtkMoveVOI()
      
   def onSliderVOI_Pt2_Z(self, evt):
      """ Fonction pour modifier la position de la VOI"""
      obj = evt.GetEventObject()
      self.vtkPipe.VOI_pt2[2] = obj.GetValue()
      self.vtkPipe.vtkMoveVOI()

   def onCloseSliderVOI(self, evt):
      """ Fermeture de la boite de dialogue SliderVOI"""
      self.sliderVOIDialog.Close(True) 

   def onButtonAllSurf(self, evt):
      """ Decoupe les surfaces connectes"""
      self.vtkPipe.aneurysmConnectFilter.SetExtractionModeToAllRegions()
      if self.vtkPipe.aneurysm<>None:
         self.vtkPipe.vtkComputeGeometricSurface()
         self.vtkPipe.vtkCreateGeometricInformations()
         if self.toolsInfoGeo.IsChecked() : 
            self.vtkPipe.vtkInfoGeo.SetInput(self.vtkPipe.GeoInfos)
      self.vtkPipe.vtkRedraw()
      
   def onButtonLargeSurf(self, evt):
      """ Conserve la surface la plus large"""
      self.vtkPipe.aneurysmConnectFilter.SetExtractionModeToLargestRegion() 
      if self.vtkPipe.aneurysm<>None:
         self.vtkPipe.vtkComputeGeometricSurface()
         self.vtkPipe.vtkCreateGeometricInformations()
         if self.toolsInfoGeo.IsChecked() : 
            self.vtkPipe.vtkInfoGeo.SetInput(self.vtkPipe.GeoInfos)
      self.vtkPipe.vtkRedraw()
      
   def onButtonArea(self, evt):
      """ Evaluation de l'aire"""
      if self.vtkPipe.aneurysm<>None:
         self.vtkPipe.vtkComputeGeometricSurface()
         self.vtkPipe.vtkCreateGeometricInformations()
         if self.toolsInfoGeo.IsChecked() : 
            self.vtkPipe.vtkInfoGeo.SetInput(self.vtkPipe.GeoInfos)
            self.vtkPipe.vtkRedraw()

   #~ def onSliderSurface(self, evt):
      #~ """ Reduction du nombre de triangle and smooth surface"""
      #~ self.sliderSurfaceDialog = wx.Dialog(self, id=wx.ID_ANY, title="Surface parameters ", size=wx.Size(220, 140), style=wx.DEFAULT_DIALOG_STYLE)
      #~ #-----------------------------------------
      #~ # Slider decimate
      #~ #-----------------------------------------
      #~ text = wx.StaticText(self.sliderSurfaceDialog, id=wx.ID_ANY, label="Decimate", pos=(10, 20))
      #~ slider0 = wx.Slider(self.sliderSurfaceDialog, id=wx.ID_ANY, value=self.vtkPipe.valDecimate, minValue=0, maxValue=100, pos=(10, 30), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      #~ slider0.Bind(wx.EVT_SCROLL, self.onSliderDecimate)
      #~ wx.StaticLine(self.sliderSurfaceDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(10, 50))
      #~ #-----------------------------------------
      #~ # Slider sagittal
      #~ #-----------------------------------------
      #~ text = wx.StaticText(self.sliderSurfaceDialog, id=wx.ID_ANY, label="Smooth surface", pos=(10, 70))
      #~ slider1 = wx.Slider(self.sliderSurfaceDialog, id=wx.ID_ANY, value=self.vtkPipe.valSmoothIter, minValue=0, maxValue=100, pos=(10, 80), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      #~ slider1.Bind(wx.EVT_SCROLL, self.onSliderSmoothIter)
      #~ wx.StaticLine(self.sliderSurfaceDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(10, 100))
      #~ #-----------------------------------------
      #~ # Bouton de fermeture de la boite de dialogue
      #~ #-----------------------------------------
      #~ close = wx.Button(self.sliderSurfaceDialog, label='Close', pos=(130, 120), size=(70, 20))
      #~ close.Bind(wx.EVT_BUTTON, self.onCloseSliderSurface)
      #~ self.sliderSurfaceDialog.Show()

   def onCloseSliderSurface(self, evt):
      """ Fermeture de la boite de dialogue SliderSurface"""
      self.sliderSurfaceDialog.Close(True)
      
   def CloseDicom2STL(self):
      """ Fermeture de l'application"""
      self.dcmUtil.remove_Images()
      self.gmshUtil.remove_Mesh()
      self.Close()
      
   def PleaseCompute(self, evt):
      """ Affichage du calcul en cour """
      evt.Skip()
      message = "Please wait, compute in progress..."
      self.busy = PBI.PyBusyInfo(message, parent=None, title="Really Busy")
      wx.Yield()

   def DestroyCompute(self):
      """ Destruction de la fenetre d'attente de fin de calcul """
      del self.busy


   def onButtonCenterline(self, evt):
      """Affiche une boite de dialogue pour le calcul des centerlines """
      
      #~ if not(self.toolsCenterline.IsChecked()):
         #~ self.vtkPipe.vtkPlotBoundingVOI()
         #~ self.meshMenu.Check(self.toolsCenterline.GetId(), True)
      
      self.centerline.listBoundary = []
      
      self.centerlineDialog = wx.Dialog(self, id=wx.ID_ANY, title="Centerlines parameters ", size=wx.Size(240, 600), style=wx.DEFAULT_DIALOG_STYLE)
      #-----------------------------------------
      # Choix de la surface support de calcul
      #-----------------------------------------
      sizer = wx.BoxSizer(wx.VERTICAL)
      self.surfaceBox = wx.RadioBox(self.centerlineDialog, id=wx.ID_ANY, label = "Choose surface", pos = (10, 20), size=wx.Size(220, 50), choices= ['Surface', 'Mesh surface'], majorDimension= 2, style = wx.RA_SPECIFY_COLS )
      #-----------------------------------------
      # Bouton de choix de la surface de calcul
      #-----------------------------------------
      selectSurf = wx.Button(self.centerlineDialog, label='Select surface', pos=(130, 80), size=(100, 25))
      selectSurf.Bind(wx.EVT_BUTTON, self.onSelectSurfCenterlineDialog)
      #-----------------------------------------
      # Choix de la surface support de calcul
      #-----------------------------------------
      wx.StaticText(self.centerlineDialog, id=wx.ID_ANY, label = "Select two boundary :", pos = (10, 120))
      self.listBoxArtery = wx.ListBox(self.centerlineDialog, id=wx.ID_ANY, pos = (20, 140), size=wx.Size(200, 160), choices= self.centerline.listBoundary, style = wx.LB_EXTENDED)
      self.listBoxArtery.Bind(wx.EVT_LISTBOX, self.onSelectContour, self.listBoxArtery)
      #~ self.listBoxArtery.Bind(wx.EVT_RIGHT_UP, self.onDeleteContour)
      self.listBoxArtery.SetSelection(0)
      #-----------------------------------------
      # Slider pour le reglage de la resolution
      #-----------------------------------------
      text = wx.StaticText(self.centerlineDialog, id=wx.ID_ANY, label="Centerline resolution", pos=(20, 320))
      slider1 = wx.Slider(self.centerlineDialog, id=wx.ID_ANY, value=self.centerline.resolution, minValue=20, maxValue=200, pos=(20, 340), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider1.Bind(wx.EVT_SCROLL, self.onSliderCenterlineResolution)
      #-----------------------------------------
      # Slider pour le reglage du filtre
      #-----------------------------------------
      text = wx.StaticText(self.centerlineDialog, id=wx.ID_ANY, label="Centerline filter", pos=(20, 380))
      slider2 = wx.Slider(self.centerlineDialog, id=wx.ID_ANY, value=self.centerline.filter, minValue=1, maxValue=10, pos=(20, 400), size=(200, -1), style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
      slider2.Bind(wx.EVT_SCROLL, self.onSliderCenterlineFilter)
      wx.StaticLine(self.centerlineDialog, id=wx.ID_ANY, style=wx.LI_HORIZONTAL, size=(200, -1), pos=(20, 440))
      #-----------------------------------------
      # Bouton de choix des frontiere de calcul
      #-----------------------------------------
      self.computeCenterLine = wx.Button(self.centerlineDialog, label='Compute CenterLine', pos=(130, 500), size=(100, 25))
      self.computeCenterLine.Bind(wx.EVT_BUTTON, self.onComputeCenterLine)
      #-----------------------------------------
      # Bouton de fermeture de la boite de dialogue
      #-----------------------------------------
      close = wx.Button(self.centerlineDialog, label='Close', pos=(130, 540), size=(100, 25))
      close.Bind(wx.EVT_BUTTON, self.onCloseCenterlineDialog)
      self.centerlineDialog.Show()

   def onSliderCenterlineResolution(self, evt):
      "User a choisi la resolution de calcul de la centerline"
      obj = evt.GetEventObject()
      self.centerline.resolution = obj.GetValue()

   def onSliderCenterlineFilter(self, evt):
      "User a choisi la resolution de calcul de la centerline"
      obj = evt.GetEventObject()
      self.centerline.filter = obj.GetValue()

   def onSelectContour(self, evt):
      "User a selectionner une arete des contours"
      namebd = self.listBoxArtery.GetSelections()
      nodeList = []
      for bd in namebd:
         for node in self.centerline.boundarys[self.centerline.listBoundary[bd]]:
            nodeList.append(node)
         
      self.centerline.createVTKContour( self.centerline.meshNode[nodeList,:])
      self.vtkPipe.vtkPlotContour(self.centerline.profileContourTubes)

   def onComputeCenterLine(self, evt):
      "User a choisi les frontieres de la centerline"
      namebd = self.listBoxArtery.GetSelections()
      
      self.centerline.computeSol(self.centerline.listBoundary[namebd[0]], self.centerline.listBoundary[namebd[1]])
      self.centerline.computeCenterLine()
      self.centerline.createVTKCenterLine('moy')
      self.vtkPipe.vtkPlotCenterline(self.centerline.profileTubes)

   def onSelectSurfCenterlineDialog(self, evt):
      "User a choisi la surface support"

      if self.surfaceBox.GetStringSelection()=="Surface":
         self.centerline.setMesh(self.vtkPipe.surfNode, self.vtkPipe.surfTable)
      if self.surfaceBox.GetStringSelection()=="Mesh surface":
         self.centerline.setMesh(self.vtkPipe.gmeshNode, self.vtkPipe.gmeshTable)
      self.centerline.getBoundary()
      self.listBoxArtery.InsertItems(self.centerline.listBoundary, 0)
      self.centerline.initializeSystem()
      self.centerline.assembleK()

   def onCloseCenterlineDialog(self, evt):
      """ Fermeture de la boite de dialogue SliderVOI"""
      self.vtkPipe.vtkDeleteCenterLine()
      self.centerlineDialog.Close(True)
      

#----------------------------------------------------------------------------------------------
class MySplashScreen(wx.SplashScreen):
   def __init__(self):
      bmp = wx.Image(os.path.join(os.path.join(os.getcwd(),'ico'),"splash.png")).ConvertToBitmap()
      wx.SplashScreen.__init__(self, bmp,
                              wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT,
                              2100, None, -1)
      self.Bind(wx.EVT_CLOSE, self.OnClose)
      self.fc = wx.FutureCall(2000, self.ShowMain)

   def OnClose(self, evt):
      # Make sure the default handler runs too so this window gets
      evt.Skip()
      self.Hide()

   def ShowMain(self):
      frame = dcm2stl(None)
      frame.Show()

#----------------------------------------------------------------------------------------------
if __name__ == "__main__":
   app = wx.App(0)
   splash = MySplashScreen()
   splash.Show()
   app.MainLoop()
   
