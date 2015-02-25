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

import string, os, logging
from logging.handlers import RotatingFileHandler

class config():
   """ Classe de fonction pour la gestion de la configuration User"""
   
   def __init__(self, path):
      " Initialisation de la classe "
      self.path = path
      self.config = {}
      self.config['TMP_PATH'] = os.path.join(os.getcwd(),'tmp')
      self.config['GMSH_PATH'] = os.path.join(os.getcwd(),'gmsh')
      self.version = 'V-1.0'
      
   def importConfig(self):
      " Lecture du fichier de configuration "
      workingPath = os.getcwd()
      os.chdir(self.path)
      file_config = open('config.txt','r')
      data = file_config.readlines()
      file_config.close()
      os.chdir(workingPath)
      for lines in data:
         line = string.split(lines)
         self.config[line[0]] = line[1]
      
   def initConfig(self):
      " Suppression de la configuration "
      self.config = {}
      self.config['TMP_PATH'] = os.path.join(os.getcwd(),'tmp')
      self.config['GMSH_PATH'] = os.path.join(os.getcwd(),'gmsh')
      
   def saveConfig(self, config):
      " Sauvegarde de la configuration "
      config.pop('TMP_PATH')
      config.pop('GMSH_PATH')
      workingPath = os.getcwd()
      os.chdir(self.path)
      file_config = open('config.txt','w')
      for data in config :
         file_config.write('%s\t%s' %(data,config[data]) + os.linesep)
      file_config.close()
      os.chdir(workingPath)
      
class logData():
   """ Classe de l'outil de log"""
   def __init__(self, path):
      " Initialisation du logger "
      
      # Repertoire de sauvegarde des logs
      self.path = path
      # création de l'objet logger qui va nous servir à écrire dans les logs
      self.logger = logging.getLogger()
      # on met le niveau du logger à DEBUG, comme ça il écrit tout
      self.logger.setLevel(logging.DEBUG)
       
      # création d'un formateur qui va ajouter le temps, le niveau
      # de chaque message quand on écrira un message dans le log
      formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
      # création d'un handler qui va rediriger une écriture du log vers
      # un fichier en mode 'append', avec 1 backup et une taille max de 1Mo
      file_handler = RotatingFileHandler(os.path.join(self.path,'activity.log'), 'a', 1000000, 1)
      # on lui met le niveau sur DEBUG, on lui dit qu'il doit utiliser le formateur
      # créé précédement et on ajoute ce handler au logger
      file_handler.setLevel(logging.DEBUG)
      file_handler.setFormatter(formatter)
      self.logger.addHandler(file_handler)
      # création d'un second handler qui va rediriger chaque écriture de log
      # sur la console
      steam_handler = logging.StreamHandler()
      steam_handler.setLevel(logging.DEBUG)
      self.logger.addHandler(steam_handler)
      
   def logMessInfo(self, msg):
      "Ecriture d'un message de log info"
      self.logger.info(msg)
      
   def logMessWarm(self, msg):
      " Ecriture d'un message error"
      self.logger.warn(msg)
   
   def logMessDebug(self, msg):
      " Ecriture d'un message error"
      self.logger.debug(msg)
   
   def logMessError(self, msg):
      " Ecriture d'un message error"
      self.logger.error(msg)
   
   def logMessCritic(self, msg):
      " Ecriture d'un message error"
      self.logger.critical(msg)   
