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
import dicom, array
import numpy as np
import os
#----------------------------------------------------------------------------------------------  
class dicomUtility():
   """ Classe de fonction pour utiliser lire les images dicom"""
   def __init__(self, write_path, log):
      """ Initialisation de la classe """
      self.write_path = write_path
      self.log = log
      
   def import_dicom(self, path, file):
      """ Fonction d'importation d'une image dicom"""
      data_dicom = dicom.read_file(os.path.join(path,file))
      return data_dicom
      
   def createPatientInformations(self, data):
      """ Creation du texte d'information patient"""
      self.PatientInfos = "--------------------------------" + os.linesep
      self.PatientInfos += "Institution name : %s" %self.get_institution_name(data) + os.linesep
      self.PatientInfos += "--------------------------------" + os.linesep
      self.PatientInfos += "Patient name : %s" %self.get_patient_name(data) + os.linesep
      self.PatientInfos += "Patient age : %s" %self.get_patient_age(data) + os.linesep
      self.PatientInfos += "Patient sex : %s" %self.get_patient_sex(data) + os.linesep
      self.PatientInfos += "Patient size : %s" %self.get_patient_size(data) + os.linesep
      self.PatientInfos += "Patient weight : %s" %self.get_patient_weight(data) + os.linesep
      self.PatientInfos += "--------------------------------" + os.linesep
      self.PatientInfos += "Body part : %s" %self.get_body_part(data) + os.linesep
      self.PatientInfos += "Sequence description : %s %s" %(self.get_modality_name(data), self.get_series_description_name(data)) + os.linesep
      self.PatientInfos += "Sequence name : %s" %self.get_sequence_name(data) + os.linesep
      self.PatientInfos += "--------------------------------" + os.linesep
      self.log.logMessInfo(self.PatientInfos)

   def get_pixel_spacing(self, data):
      """Recupere la distance entre pixel"""
      a = data[0x28,0x30].value
      return float(a[0]),float(a[1])
      
   def get_pixel_codage(self, data):
      """Recupere le codage de l'image"""
      return data[0x28,0x101].value
   
   def get_slice_spacing(self, data):
      """Recupere l'espacement entre les coupes"""
      return float(data[0x18,0x50].value)

   def get_slice_thickness(self, data):
      """Recupere l'espacement entre les coupes"""
      return float(data[0x18,0x50].value)

   def get_slice_position(self, data):
      """Recupere la position d'une coupe"""
      a = data[0x20,0x32].value
      return np.array([float(a[0]),float(a[1]),float(a[2])])
   
   def get_slice_number(self, data):
      """Recupere l'indice de la coupe"""
      return int(data[0x20,0x13].value)
   
   def get_slice_data(self, data):
      """Recupere les donnees d'une coupe"""
      try :
         val = data.pixel_array
      except :
         val = np.zeros((data.Rows,data.Columns),int)
         self.log.logMessError("DICOM PROBLEM : pixel_array does not exist - dicom file is compressed.")
      #~ values = np.fromstring(data.PixelData,dtype=np.int16)
      #~ values = values.reshape((data.Rows,data.Columns))
      return val
      
   def get_slice_size(self, data):
      """Recupere la resolution d'une coupe"""
      return data[0x28,0x10].value,data[0x28,0x11].value
   
   def get_max_pixel_value(self, data):
      """Recupere la valeur maxi d'une coupe"""
      try :
         value = data[0x28,0x107].value
      except :
         value = 2**data[0x28,0x101].value
      return value
   
   def get_min_pixel_value(self, data):
      """Recupere la valeur mini d'une coupe"""
      try :
         value = data[0x28,0x106].value
      except :
         value = 0
      return value
   
   def get_patient_name(self, data):
      """Recupere les informations """
      try :
         value = data[0x10,0x10].value
      except :
         value = 'not defined'
      return value

   def get_patient_birthdate(self, data):
      """Recupere les informations """
      try :
         value = data[0x10,0x30].value
      except :
         value = 'not defined'
      return value

   def get_patient_sex(self, data):
      """Recupere les informations """
      try :
         value = data[0x10,0x40].value
      except :
         value = 'not defined'
      return value

   def get_patient_age(self, data):
      """Recupere les informations """
      try :
         value = data[0x10,0x1010].value
      except :
         value = 'not defined'
      return value

   def get_patient_size(self, data):
      """Recupere les informations """
      try :
         value = data[0x10,0x1020].value
      except :
         value = 'not defined'
      return value

   def get_patient_weight(self, data):
      """Recupere les informations """
      try :
         value = data[0x10,0x1030].value
      except :
         value = 'not defined'
      return value
      
   def get_body_part(self, data):
      """Recupere les informations """
      try :
         value = data[0x18,0x15].value
      except :
         value = 'not defined'
      return value

   def get_institution_name(self, data):
      """Recupere les informations """
      try :
         value = data[0x8,0x80].value
      except :
         value = 'not defined'
      return value

   def get_modality_name(self, data):
      """Recupere les informations """
      try :
         value = data[0x8,0x60].value
      except :
         value = 'not defined'
      return value

   def get_sequence_name(self, data):
      """Recupere les informations """
      try :
         value = data[0x18,0x24].value
      except :
         value = 'not defined'
      return value

   def get_series_description_name(self, data):
      """Recupere les informations """
      try :
         value = data[0x8,0x103e].value
      except :
         value = 'not defined'
      return value

   def add_slice2vtk(self, data, k) :
      """add a slice to vti volum with vtk"""
      data.astype('int16').tofile(os.path.join(self.write_path,'image.raw.' + str(k)))
   
   def remove_Images(self):
      """ Suppression des images de travail """
      images = os.listdir(self.write_path)
      for image in images:
         try : 
            os.remove(os.path.join(self.write_path,image))
         except :
            self.log.logMessInfo("DICOM PROBLEM : Inexist file : "+image)

