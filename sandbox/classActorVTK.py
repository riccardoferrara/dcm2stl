# -*- coding: utf-8 -*-
#!===========================================================================
#!
#! Copyright 2000-2015 CNRS.
#!
#! This file is part of a software (dicom2stl) which is a computer program 
#! which purpose is to modelize interaction problems (contact, multi-Physics,etc).
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

import vtk,wx

class actorPolyDataVTK():
   """ Classe d'objet vtk pour rendu 3D """
   
   def __init__(self, ren, widget):
      " Initialisation de la classe "
      self.ren    = ren
      self.widget = widget
      self.vtkActor = None
      self.color = (1, 1, .9412)
      self.opacity = 1.0
      
   def add(self):
      "fonction ajout du rendu"
      self.ren.AddActor(self.vtkActor)Mapper
      self.widget.Render()
      
   def redraw(self):
      "fonction de rendu"
      self.widget.Render()
      
   def delete(self):
      "fonction de suppression"
      self.ren.RemoveActor(self.vtkActor)
      self.widget.Render()
      self.centerlineActor = None
      
   def createPolyData(self, polydata):
      "fonction de creation de la visualisation"
      self.vtkMapper = vtk.vtkPolyDataMapper()
      if vtk.VTK_MAJOR_VERSION <= 5:
          self.vtkMapper.SetInput(polydata.GetOutput())
      else:
          self.vtkMapper.SetInputConnection(polydata.GetOutputPort())
      
      self.vtkActor = vtk.vtkActor()
      self.vtkActor.SetMapper(self.vtkMapper)
      self.vtkActor.GetProperty().SetDiffuseColor(self.color[0], self.color[1], self.color[2])
      self.vtkActor.GetProperty().SetSpecularColor(1, 1, 1)
      self.vtkActor.GetProperty().SetSpecular(.4)
      self.vtkActor.GetProperty().SetSpecularPower(50)
      self.vtkActor.GetProperty().SetOpacity(self.opacity)
      
   def setColor(self):
      "fonction de choix de la couleur"
      dialog = wx.ColourDialog(self)
      dialog.GetColourData().SetChooseFull(True)
      
      if dialog.ShowModal() == wx.ID_OK:
         data = dialog.GetColourData()
         color = data.GetColour().Get()
         self.color = (color[0]/255.0, color[1]/255.0, color[1]/255.0)
         self.vtkActor.GetProperty().SetDiffuseColor(self.color[0], self.color[1], self.color[2])
         self.redraw()
         
      dialog.Destroy()
