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

import os
import string
import time
import subprocess
from sys import platform as _platform

class gmshUtility():
   """ Fonction pour utiliser les routines gmsh de l'exterieur """
   
   def __init__(self, path, log):
      """ Initialisation de la classe """
      self.path = path
      self.log = log
      
      self.remesh  = ' -2 format vtk -o SurfaceMesh.vtk Remeshing.geo > remesh.log'
      self.extrude = ' -3 format vtk -o SolidMesh.vtk Extrude.geo > extrude.log'
      
      if _platform == "linux" or _platform == "linux2":
          # linux
          self.process = './gmsh'
      elif _platform == "darwin":
          # OS X
          self.process = 'gmsh'
      elif _platform == "win32":
          # Windows...
          self.process = 'gmsh.exe'
      
   def Remeshing(self):
      """ Remeshing with gmsh """
      workingPath = os.getcwd()
      os.chdir(self.path)
      self.log.logMessInfo("GMSH COMMAND : " + self.process + self.remesh)
      p = subprocess.Popen(self.process + self.remesh,stdin=None,stdout=None, shell=True)
      p.wait()
      os.chdir(workingPath)
      
   def remove_Mesh(self):
      """ Suppression des fichiers temporaires du maillage """
      try : 
         os.remove(os.path.join(self.path, 'Surface.stl'))
      except :
         self.log.logMessInfo("GMSH PROBLEM : Inexist file : Surface.stl")
      
      try : 
         os.remove(os.path.join(self.path, 'SurfaceMesh.vtk'))
      except :
         self.log.logMessInfo("GMSH PROBLEM : Inexist file : SurfaceMesh.vtk")
      
      try : 
         os.remove(os.path.join(self.path, 'SolidMesh.vtk'))
      except :
         self.log.logMessInfo("GMSH PROBLEM : Inexist file : SolidMesh.vtk")
         
   def Extrude(self):
      """ Extrude mesh from  surface """
      workingPath = os.getcwd()
      os.chdir(self.path)
      self.log.logMessInfo("GMSH COMMAND : " + self.process + self.extrude)
      p = subprocess.Popen(self.process + self.extrude,stdin=None,stdout=None, shell=True)
      p.wait()
      os.chdir(workingPath)
      
