#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: jling
Mayerlab bids utils
"""

import os.path, glob
from mayerlab.preprocessing import utils
    
    
class Bidsmake:
    """
    Simplify converted dicom into BIDS dirctory structure

    All arguments are required but dicom arguments setup to either accept a known path 'dicom_full_path'
    or determine path from a root path and a search string. In search case, last found item matching
    the string is used under assumption that the last acquistion is the good one.
    Intialization of object will report dicom path and full test string of what it would want crete.

    'modality' conforms to BIDS specification [anat,func,dwi]
    'label' is string identifing acquistion. This will be built out if 'modality=func'
    'session_id' = how the session is labeled in bids tree
    'raw_path_root' = bids tree path; should usually be /xx/xx/xx/rawdata
    
    """

    
    
    def __init__(self, subject_id=None, session_id=None, modality=None, label=None, raw_path_root=None, dicom_path_root=None, dicom_search_str=None, dicom_full_path=None, json_fields=None):
        
        if subject_id == None  or modality == None or label == None or raw_path_root == None:
            raise ValueError("BIDSMAKE: subject_id, modality, label, and raw_path_root are required.")

        # check for valid modality
        if modality not in ['anat', 'func', 'dwi', 'fmap', 'perf', 'meg', 'eeg', 'ieeg', 'beh']:
            raise ValueError("BIDSMAKE: " + modality + " not valid modality(Data type)")

        self.subject_id = "sub-"+subject_id
        self.o_file_name = self.subject_id

        if session_id != None:
            if isinstance(session_id, int):
                #self.session_id = "ses-"+str(session_id).zfill(3)
                self.session_id = "ses-"+str(session_id)
            else:
                self.session_id = "ses-"+session_id
                
            self.o_file_name = self.o_file_name + "_" + self.session_id
                
        else:
            self.session_id = ""
            

        
        # label
        self.modality = modality
        if modality == "func":
            self.label = "task-" + label
        elif modality == "dwi":
            self.label = label + "_dwi"
        elif modality == "fmap":
            self.label = label + "_epi"
        else:
            # checks bids-apps like MRIQC, underscore characters in the key value broke MRIQC processing
            # what does the BIDS spec say about this??????
            self.label = label
        
        
        self.o_file_name = self.o_file_name + "_" + self.label

        
        if dicom_full_path is None:
            dicom_input_list = glob.glob(os.path.join(dicom_path_root, dicom_search_str))
            if len(dicom_input_list):
            
                # glob order can be arbitrary depending on how file objects were created, 
                dicom_input_list = sorted(dicom_input_list)
                # use the last item from the sorted list as our dicom input path                
                self.dicom_full_path = dicom_input_list[-1]
                
            else:
                raise ValueError("Could not determine dicom source path:\n  dicom_path_root:", dicom_path_root, "\n  dicom_search_str:", dicom_search_str, "\n  dicom_full_path:", dicom_full_path)
        else:
            self.dicom_full_path = dicom_full_path
            
            
        
        self.dicom_path_root = dicom_path_root
        self.raw_path_root = raw_path_root
        self.dicom_search_str = dicom_search_str
        
        self.bids_out_path = os.path.join(raw_path_root, self.subject_id, self.session_id, self.modality)
        self.o_file_test_path = os.path.join(self.bids_out_path, self.o_file_name + ".nii.gz")
        
        self.json_fields = json_fields
        
        
    
        
    def __repr__(self):
        return "BIDSMAKE DICOM SRC:  " + self.dicom_full_path + "\nBIDSMAKE TARGET: " + self.o_file_test_path
        

    def dicom_to_nifti(self):
        """
        usally runs conversion and alters storage
        """        
        if not os.path.isfile(self.o_file_test_path) and self.dicom_full_path is not None:
            utils.prettyOut("BIDSMAKE: Convert Dicom to Nifti...")
            print(self)
            os.makedirs(self.bids_out_path, exist_ok=True)
            utils.Cmd("dcm2niix -l o -b y -d 1 -z y=pigz -o "+ self.bids_out_path +" -f "+ self.o_file_name +" " + self.dicom_full_path).run()
            
            
            if len(self.json_fields):                
                import json
                import glob
                
                # find json file from conversions
                json_file_test_path = os.path.join(self.bids_out_path, self.o_file_name + ".json")
                
                # add entry in json file and save out
                with open(json_file_test_path, 'r+') as f:
                    data = json.load(f)

                    # for every json element in list add a json element
                    for json_field in self.json_fields:
                        for i in json_field:
                            # if the value is a string and contains "*" expand the path to the relative session path
                            if isinstance(json_field[i], str) and ('*' in json_field[i]):
                                intend_file = glob.glob(os.path.join(self.bids_out_path, json_field[i]))
                                if len(intend_file):
                                    data[i] = "/".join(intend_file[-1].split("/")[-3:])

                            else:
                                data[i] = json_field[i]

                    f.seek(0)
                    json.dump(data, f, indent=4)
                    f.truncate() 
  
        else:
            utils.prettyOut("BIDSMAKE: Nifti output already exists")
            print(self)
    



    @classmethod
    def get_bids_name(self, subject_id, session_id=None, label=None):
        """
        returns a bids filename from elements
        """
        file_name = "sub-"+subject_id
        
        if session_id:
            file_name = file_name+"_ses-"+session_id
            
        if label:
            file_name = file_name+"_"+label
            
        return file_name+".nii.gz"
    
    
    
    
    @classmethod
    def create_description_file(self, raw_path_root=None, study_name=None, author_list=[]):
        """
        Writes the required description file to raw path
        """
        import json
        
        if raw_path_root == None  or study_name == None or not len(author_list):
            raise ValueError("BIDSMAKE.create_description_file(): raw_path_root, study_name, author_list are required.")        
        
        
        if not os.path.isdir(raw_path_root):
            os.makedirs(raw_path_root, exist_ok=True)
            
            
        for_json = {
          "Name": study_name,
          "BIDSVersion": "1.2.2",
          "License": "CC0",
          "DatasetType":"raw", 
          "Authors": author_list
        }
        
        with open( os.path.join(raw_path_root, 'dataset_description.json'), "w" ) as write_file:
            json.dump(for_json, write_file)


        
#    @classmethod
#    def get_bids_path(self, ):
#        """
#        returns a bids filepath from elements
#        """

class SideCar(object):

    
    """
    load some value for json sidecar key
    """
    @staticmethod
    def get_param(json_path=None, param_key=None):
        
        if json_path is not None:
            if not os.path.isfile(json_path):
                raise ValueError("Sidecar path not found")
        
        
            import json

            with open(json_path, 'r') as json_file:
                data = json.load(json_file)
                return (data[param_key])
                        
                
            
    
    @staticmethod
    def set_param(json_path=None, param_key=None, param_value=None):
        
        if json_path is not None:
            if not os.path.isfile(json_path):
                raise ValueError("Sidecar path not found")
        
        
            import json

            with open(json_path, 'r') as json_file:
                data = json.load(json_file)
                data[param_key] = param_value
                json_file.seek(0)
                json.dump(data, json_file, indent=4)
                json_file.truncate() 
                
            
            
                    

            