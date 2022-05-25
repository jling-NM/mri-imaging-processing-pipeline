# -*- coding: utf-8 -*-

"""@package utils

Created on Wed Dec  2 13:28:49 2015

@author: josef ling

"""

import sys
import errno
import os

class Ursi:

    def __init__(self, subject_id, id_prefix="M871"):

        if len(subject_id) == 9 :
            self.full     = subject_id
            self.short     = subject_id[4:9]
        elif len(subject_id) == 5 :
            self.full     = id_prefix + subject_id
            self.short     = subject_id
        else :
            self.full     = subject_id
            self.short     = subject_id




class Study:

    def __init__(self, study_label):

        self.label           = study_label.upper()
        self.settings        = self.get_settings(self.label)
        if self.settings == None:
            raise ValueError("Study label '" + self.label + "' not found in available studies.")
            
        self.path            = self.settings['dataPath']
        self.bidsPath        = self.settings['bidsPath']
        self.freesurferPath  = self.settings['freesurferDataPath']
        self.mriqcPath       = self.settings['mriqcPath']
        
        
    @property
    def protocol(self):
        return self.settings['protocol']

    @property
    def dataPath(self):
        return self.settings['dataPath']
    
    @property
    def archivePath(self):
        return self.settings['archivePath']
    
    @property
    def paramsPath(self):
        return self.settings['paramsPath']
    
    
    def get_settings(self, label):

        import sys, json
        
        try:
            with open('/export/research/analysis/human/amayer/shared/apps/python/mayerlab/preprocessing/study_settings.json') as json_file:
                data = json.load(json_file) 
                
                for study in data['study']:
                    if study['label'] == label:
                        return study
                    
            return None
        
        except IOError as e:
            print(("I/O error({0}): {1}".format(e.errno, e.strerror)))
            raise
        except:
            print(("Unexpected error:", sys.exc_info()[0]))
            raise







    #    takes full run type name and returns data specific to the run number included in that name
    def get_run_label_settings(self,run_label):

        run_settings_dict = dict(run_count=None,time_images=None,num_slices=None,TR=None)

        import sys

        try:

            import re

            run_type = run_label.rstrip("0123456789")

            #    if no number in run_label we'll assume run1 settings data
            run_num  = 666
            if len(list(map(int, re.findall(r'\d+', run_label)))):
                run_num  = list(map(int, re.findall(r'\d+', run_label)))[-1]


            settings_file_path = self.get_settings_path()
            settings_file = open(settings_file_path, "r")

            try:
                for line in settings_file:

                    if line.startswith( "set " + run_type + "_run_count" ):
                        d_list = line.split("=")[1].strip().translate(str.maketrans('', '', '()[]')).split()
                        if run_num > len(d_list):
                            run_settings_dict['run_count'] = d_list[0]
                        else:
                            run_settings_dict['run_count'] = d_list[run_num-1]

                        #if(int(run_num) > int(run_settings_dict['run_count']) ):
                        #    raise ValueError("Requested run number exceeds run count")


                    if line.startswith( "set " + run_type + "_time_images" ):
                        d_list = line.split("=")[1].strip().translate(str.maketrans('', '', '()[]')).split()
                        if run_num > len(d_list):
                            run_settings_dict['time_images'] = d_list[0]
                        else:
                            run_settings_dict['time_images'] = d_list[run_num-1]


                    if line.startswith( "set " + run_type + "_num_slices" ):
                        d_list = line.split("=")[1].strip().translate(str.maketrans('', '', '()[]')).split()
                        if run_num > len(d_list):
                            run_settings_dict['num_slices'] = d_list[0]
                        else:
                            run_settings_dict['num_slices'] = d_list[run_num-1]

                    if line.startswith( "set " + run_type + "_TR" ):
                        d_list = line.split("=")[1].strip().translate(str.maketrans('', '', '()[]')).split()
                        if run_num > len(d_list):
                            run_settings_dict['TR'] = d_list[0]
                        else:
                            run_settings_dict['TR'] = d_list[run_num-1]

            finally:
                settings_file.close()

            return run_settings_dict


        except ValueError as e:
            print(("Input error({0}): {1}".format(run_label, e)))
            raise
        except IOError as e:
            print(("I/O error({0}): {1}".format(e.errno, e.strerror)))
            raise
        except:
            print(("Unexpected error:", sys.exc_info()[0]))
            raise



    def create_bids_description_file(self):
        from mayerlab.preprocessing import bids
        
        bids.Bidsmake.create_description_file(raw_path_root=self.bidsPath, study_name=self.label, author_list=['Andrew Mayer'])
        



    def __str__(self):
        """
        override/define print response
        """
        return "\nSTUDY\n\t LABEL:"+self.label+"\n\t PATH:" + self.path + "\n"






class Subj:

    def __init__(self, subject_id, study, visit=1, id_prefix="M871", dx_int=None):

        self.ursi             = Ursi(subject_id, id_prefix)
        self._study           = study
        self._dx_int          = dx_int

        if visit == None:
            self._visit       = None
        else:
            self._visit       = visit






    @property
    def bids_path(self):
        return os.path.join(self.study.bidsPath, self.visit_path)

    @property
    def bids_id(self):
        return 'sub-'+self.ursi.short
    
    @property
    def bids_visit(self):
        return 'ses-'+ str(self.visit)
    
    def bids_modality_path(self, bidsType=''):
        return os.path.join(self.study.bidsPath, self.visit_path, bidsType)
    
    def derivatives_path(self, bidsType=None):
        return os.path.join(self.study.dataPath, self.visit_path, bidsType)
    
    
    def bids_file_path(self, bidsType=None, bidsLabel=None, task=None, run=None, extension='nii.gz'):
        
        if bidsType == 'func':
            filename_suffix = 'epi'

        elif bidsType == 'fmap':
            filename_suffix = 'epi'

        elif bidsType == 'dwi':
            filename_suffix = bidsType
            
        elif bidsType == 'perf':
            filename_suffix = ""
            
        else:
            filename_suffix = ""
            

        if task != None:
            
            if run == None:
                run = 1
            
            file_path = os.path.join(self.study.bidsPath, 
                                     self.visit_path,
                                     bidsType,
                                     "_".join([str(x) for x in [self.bids_id, "ses-"+str(self._visit), 'task-' + task.lower(), 'run-' + str(run), bidsLabel + "." + extension] if len(x)])
                                     )

        else:
            
            file_path = os.path.join(self.study.bidsPath, 
                                      self.visit_path,
                                      bidsType,
                                      "_".join([str(x) for x in [self.bids_id, "ses-"+str(self._visit), filename_suffix, bidsLabel + "." + extension] if len(x)])
                                      )

            

        if not os.path.isfile(file_path):
            file_path = None
        
        return file_path




    def bids_file_join(self, *parts):
        s = "_".join([self.bids_id, self.bids_visit])
        for part in parts:
            s = "_".join([s, part])

        return s



    
    @property
    def visit(self):
        return self._visit

    # @property
    # def visit_num(self):
    #     num = [int(s) for s in list(self._visit) if s.isdigit()][0]
    #     return num

    @property
    def visit_path(self):
        if self._visit == None:
            return "sub-"+self.ursi.short
        else:
            return "sub-"+self.ursi.short+"/"+"ses-"+str(self._visit)
    
    @property
    def visit_settings_path(self):    
        return "/".join([self.study.dataPath, self.visit_path, 'visit_description.json'])
    
    
    @visit.setter
    def visit(self, value):
        self._visit = value

    @property
    def dx_int(self):
        return self._dx_int

    @dx_int.setter
    def dx_int(self, value):
        self._dx_int = value

    @property
    def study(self):
        return self._study

    @property
    def sessions(self):
        """
        Get list of all sessions for this subject in this study
        
        """
        
        import glob
        import os.path as op
        
        sessions = glob.glob(op.join(self.study.archivePath, self.ursi.full, 'Study*'))

        return sorted(sessions)


    def get_run_path_list(self, sessions, protocol, run_num):
        """
        Get list of all runs for this protocol acqusition
        
        """
        
        import glob
        import os.path as op

        runs = []
        for session in sessions:
            if protocol['bidsType'] == 'func':
                if protocol['bidsLabel'] == 'sbref':
                    session_runs = glob.glob(op.join(self.study.archivePath, self.ursi.full, session, protocol['protocolName'] + '_run'+str(run_num)+'_AP_SBRef_0???'))
                else:
                    session_runs = glob.glob(op.join(self.study.archivePath, self.ursi.full, session, protocol['protocolName'] + '_run'+str(run_num)+'_AP_0???'))

            else:
                session_runs = glob.glob(op.join(self.study.archivePath, self.ursi.full, session, protocol['protocolName'] + '_0???'))

            # add what is found in each session
            runs.extend(session_runs)

        return sorted(runs)
    
    
    
    
    def visit_settings_to_bids(self, subj_visit_settings):
        """
        Use visit mapping to init bids rawdata for this subject
        """
        from mayerlab.preprocessing import bids
        
        for run_protocol in subj_visit_settings["protocol"]:
            
            if not len(run_protocol['srcDirectory']):
                print("\n+++ No data found for " + run_protocol['label'])
            else:
                
                # now convert to bids format
                if run_protocol['bidsType'] == 'anat':
                    run_label = run_protocol['bidsLabel']

                elif run_protocol['bidsType'] == 'fmap':
                    run_label = '_'.join(['acq-'+run_protocol['label'].lower(), run_protocol['bidsLabel'] ])
                    
                elif run_protocol['bidsType'] == 'dwi':
                    run_label = '_'.join([ run_protocol['bidsLabel'], 'run-' + str(run_protocol['runNumber'])])

                elif run_protocol['bidsType'] == 'perf':
                    run_label = '_'.join([ 'acq-' + run_protocol['label'].lower(), run_protocol['bidsLabel'] ])
                    
                else:
                    run_label = '_'.join([ run_protocol['label'].lower(), 'run-' + str(run_protocol['runNumber']), run_protocol['bidsLabel'] ])
                    
                    
                
                # convert it
                bids.Bidsmake(subject_id=self.ursi.short, 
                    session_id=self.visit, 
                    modality=run_protocol['bidsType'], 
                    label=run_label,
                    raw_path_root=self.study.bidsPath, 
                    dicom_full_path=run_protocol['srcDirectory'],
                    json_fields=run_protocol['jsonFields'] ).dicom_to_nifti()   
                        
                        
                

        
    def read_visit_settings(self):
        """
        Read subject visit json file
        """
        import os.path
        import json
        
        try:
            if os.path.isfile(self.visit_settings_path):
                print("   +++ Read subject visit settings +++")
                with open(self.visit_settings_path, 'r') as fh:
                    s = json.load(fh)
                
                return s
            
            else:
                return None
            
        except IOError as e:
            print(("I/O error({0}): {1}".format(e.errno, e.strerror)))
            raise
        except:
            print(("Unexpected error:", sys.exc_info()[0]))
            raise
        
        
        
        
        
    def write_visit_settings(self, subj_visit_settings_data_dict):
        """
        Write subject visit json file
        """
        import json
        import os
        
        try:
            os.makedirs(os.path.dirname(self.visit_settings_path), 0o770)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        
        print("   +++ Write subject visit settings +++")
        with open(self.visit_settings_path, 'w') as fh:
            json.dump(subj_visit_settings_data_dict, fh, indent=4)
                
               
                
                
    def get_subj_task_map(self):
        """
        Show me which EPI task and how many runs available for this subject.
        This uses the visit json file for the subject but could also be made
        BIDS-app compatible by reading rawdata folder instead. i.e. the rawdata
        folder becomes the test of availability.

        Returns
        -------
        proc_map : TYPE
            DESCRIPTION.

        """
        try:
            proc_map = {}
            subj_visit_settings = self.read_visit_settings()
            if subj_visit_settings == None:
                raise ValueError("Subject {} not yet intialized? No visit {} settings found.".format(self.ursi.full, self.visit))
            
            for run_protocol in subj_visit_settings["protocol"]:
                if (run_protocol['bidsType'] == 'func') and (run_protocol['bidsLabel'] == 'bold') and len((run_protocol['srcDirectory'])):
                    if run_protocol['label'].lower() in proc_map:
                        proc_map[run_protocol['label'].lower()] = int(proc_map[run_protocol['label'].lower()]) + 1
                    else:
                        proc_map[run_protocol['label'].lower()] = 1
            
            return proc_map
    
        except:
            raise
            
            
    

    def study_to_visit_map(self, disableProtocolChecks=False):
        """
        """
        
        import glob
        from mayerlab.preprocessing import gui
        from datetime import datetime
        
        # a description file is required at the rawdata root
        if not os.path.isfile( os.path.join(self.study.bidsPath, 'dataset_description.json') ):
            self.study.create_bids_description_file()
        
        
        # Do we already have init file?
        subj_visit_settings = self.read_visit_settings()

        if subj_visit_settings == None:
            
            # get sessions from archive
            if not len(self.sessions):
                raise ValueError("No sessions for subject. Check URSI")
            
            if len(self.sessions) == 1:
                visit_sessions = self.sessions
            else:
                # use visit to preselect a session for user
                preselected_item = (self.visit - 1)
                user_selected_sessions_index = gui.get_user_selections('Select Sessions', 'Select the session(s) which provide data for visit ' + str(self.visit), [s.split('/')[-1] for s in self.sessions], preselectedList=[preselected_item])
                # build list of selected sessions for dicom search
                visit_sessions = []
                for i in user_selected_sessions_index:
                    visit_sessions.append(self.sessions[i])
                
            if not len(visit_sessions):
                raise ValueError("No session selected")
                
                

            
            # this will be out visit map
            subj_visit_settings = {"dateGenerated": datetime.now().strftime("%Y%m%d %H:%M:%S"), "userGenerated": os.environ['USER'], "subjectID": self.ursi.full, "visit": self.visit, "protocol": [] }
    
            # kludgy fix for dwi to have run numbers so files can be in order
            dwi_counter = 0
            
            # go through protocol looking for each item in protocol
            for run_protocol in self.study.protocol:
                
                for run_number in range(1,len(run_protocol['imageNum'])+1):

                    visit_protocol = {}
                    visit_protocol['label'] = run_protocol['label']
                    visit_protocol['bidsType'] = run_protocol['bidsType']
                    visit_protocol['bidsLabel'] = run_protocol['bidsLabel']
                    visit_protocol['jsonFields'] = run_protocol['jsonFields']
                    #print("run num", run_number, "of", len(run_protocol['imageNum']))

                    # kludgy fix for dwi to have run numbers so files can be in order; assumes study settings has correct run order
                    if run_protocol['bidsType'] == 'dwi':
                        dwi_counter = dwi_counter + 1
                        visit_protocol['runNumber'] = dwi_counter                        
                    else:
                        visit_protocol['runNumber'] = run_number
                    
                    
                    run_paths = self.get_run_path_list( visit_sessions, run_protocol, run_number )
                    
                    if len(run_paths) == 0:
                        
                        visit_protocol['srcDirectory'] = []
                        
                        if not disableProtocolChecks:
                            msg_missing_data = "Missing " + visit_protocol['bidsType'] + " data for " + visit_protocol['label'] + " " + visit_protocol['bidsLabel'] + " run " + str(run_number) + "\n"
                            confirm_protocol_skip = gui.get_user_confirm('Confirm Missing Data', msg_missing_data + "\n\nIf you want to disable ALL study protocol checks (not recommended)\nyou can add an additional flag like this:\npreprocessing.1.init.subject.visit.py -v "+str(self.visit)+" -s "+ self.ursi.short +" --disableProtocolChecks \n\nContinue now without this data?")
                            if not confirm_protocol_skip:
                                raise ValueError("+++ " + msg_missing_data)
                        
                    elif len(run_paths) > 1:
                        
                        # use visit to preselect a session for user
                        preselected_item = len(run_paths) - 1
                        selected_index = gui.get_user_selections('Select Run', 'Select the '+run_protocol['label'] + ' ' +run_protocol['bidsType']+' run which provides data for visit ' + str(self.visit), ["/".join(s.split("/")[-2:]) for s in run_paths], preselectedList=[preselected_item])
                        
                        if not len(selected_index):
                            raise ValueError("Run selection cancelled or no selection made")

                        if len(selected_index) > 1:
                            raise ValueError("You may select only one run")

                            
                        if not disableProtocolChecks:
                            if len(list(glob.iglob(os.path.join(run_paths[selected_index[0]], '*.dcm')))) != run_protocol['imageNum'][run_number-1]:
                                raise ValueError("Unexpected number of dicom files in " + run_paths[selected_index[0]])
            
                        
                        visit_protocol['srcDirectory'] = run_paths[selected_index[0]]
                        
                    else:
                        visit_protocol['srcDirectory'] = run_paths[0]
                        
                    
                    # update map
                    subj_visit_settings["protocol"].append(visit_protocol)
    
    

            # write each out map in json format or just let bids format do the rest?
            self.write_visit_settings(subj_visit_settings)
            
            
            
        
        # With settings map create or read let's convert dicom to BIDS
        self.visit_settings_to_bids(subj_visit_settings)




        


                
    def __str__(self):
        """
        override/define print response
        """
        return "SUBJ\n URSI:"+self.ursi.full+"\n STUDY:"+self.study.label+"\n VISIT:"+str(self.visit)+"\n DX CODE:"+str(self.dx_int)+"\n\nUse pprint(vars(subj)) for more detail\n"







class BIDSStudy:
    """
    Study class uses central json file for intializing study parameters.
    
    This BIDSStudy class is a reduce class for operating inside a BIDS-app where the user inputs the path and data is in BIDS format
    
    """
    import os
    
    def __init__(self, dataPath, bidsPath, studyLabel=""):

        self.label           = studyLabel
        self.path            = dataPath
        self.bidsPath        = bidsPath
        
        
    @property
    def dataPath(self):
        return self.path


    def __str__(self):
        """
        override/define print response
        """
        return "\nSTUDY\n\t LABEL:"+self.label+"\n\t PATH:" + self.path + "\n"











class Env:

    runtime_path = '/export/research/analysis/human/amayer/shared/apps'
    
    from os import environ
    in_container_value = environ.get('MAYERLAB_IN_CONTAINER')
    if in_container_value == 'True':
        runtime_path = '/apps'
        
    
    #    static or class attribute as opposed to instance attribute
    afni_path           = runtime_path + "/afni_current_64"
    fsl_path            = runtime_path + "/fsl"
    c3d_path            = runtime_path + "/c3d"
    shared_script_path  = runtime_path + "/shared_proc_scripts"
    barebrain_path      = runtime_path + "/barebrain"
    standard_masks_path = runtime_path + "/standard_masks"
    freesurfer_path     = runtime_path + "/freesurfer"

        
    # template map
    template_lib = {}
    template_lib['mni_icbm152_nlin_asym_2009c'] = {}
    template_lib['mni_icbm152_nlin_asym_2009c']['head'] = runtime_path + '/templates/mni_icbm152_nlin_asym_2009c/mni_icbm152_t1_tal_nlin_asym_09c.nii.gz'
    template_lib['mni_icbm152_nlin_asym_2009c']['brain'] = runtime_path + '/templates/mni_icbm152_nlin_asym_2009c/mni_icbm152_t1_tal_nlin_asym_09c_brain.nii.gz'
    template_lib['mni_icbm152_nlin_asym_2009c']['mask'] = runtime_path + '/templates/mni_icbm152_nlin_asym_2009c/mni_icbm152_t1_tal_nlin_asym_09c_brain_mask.nii.gz'
    template_lib['mni_icbm152_nlin_asym_6gen'] = {}
    template_lib['mni_icbm152_nlin_asym_6gen']['head'] = runtime_path + '/fsl/data/standard/MNI152_T1_1mm.nii.gz'
    template_lib['mni_icbm152_nlin_asym_6gen']['brain'] = runtime_path + '/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz'
    template_lib['mni_icbm152_nlin_asym_6gen']['mask'] = runtime_path + '/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii.gz'
    template_lib['TT_N27'] = {}
    template_lib['TT_N27']['head'] = runtime_path + '/templates/TT_N27+tlrc.HEAD'
    template_lib['TT_N27']['brain'] = runtime_path + '/templates/TT_N27+tlrc.HEAD'
    template_lib['TT_N27']['mask'] = runtime_path + '/templates/TT_N27+tlrc.HEAD'    
    
    
    
    # how much memory install on linux system
    @staticmethod
    def inst_ram_gb():
        # parse the whole table into a dictionary
        meminfo = dict((i.split()[0].rstrip(':'),int(i.split()[1])) for i in open('/proc/meminfo').readlines())
        # return what you want from that
        return meminfo['MemTotal']/1024/1024

    @staticmethod
    def inst_ram_mb():
        # parse the whole table into a dictionary
        meminfo = dict((i.split()[0].rstrip(':'),int(i.split()[1])) for i in open('/proc/meminfo').readlines())
        # return what you want from that
        return meminfo['MemTotal']/1024
    
    @staticmethod
    def freesurfer_version():
        """
        Return freesurfer version as reported by recon-all binary
        This will allow users to map their own fs paths
        """
        import subprocess
        from os.path import sep
        from os.path import join as opj
        
        fs_version = ''
        
        subproc_res = subprocess.run('which recon-all', shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        p = sep.join(subproc_res.stdout.strip().split('/')[0:-2])
        with open(opj(p, 'build-stamp.txt')) as fh:
            fs_version = fh.readline().strip()
        
        return fs_version



    @staticmethod
    def running_in_container():
        """
        Need a way of knowing if we are running in a container.
        """
        from os import environ
        in_container_value = environ.get('MAYERLAB_IN_CONTAINER')
        
        if in_container_value == 'True':
            return True
        else:
            return False
        
        
        
        

# for raising errors when calling external programs
class CmdProcError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class Log:
    """
    """
    file_path = None
    verbose = None

    
    def __init__(self, log_path="", verbose=False):
        self.file_path = log_path
        self.verbose = verbose
        

    def write(self, txt=""):
        if self.file_path is not None:
            with open(self.file_path, 'a') as f:
                f.write(txt+"\n")


class Cmd:

    def __init__(self, cmd_str="", verbose=True, log=None):
        # what to run
        self.cmd_str = cmd_str
        # if invoking afni gui; default to same detach as in terminal
        self.afni_detach  = False
        # show stuff running or just echo the command called?
        self.verbose = verbose
        # will stdout by captured to a file
        self.log = log


    def run(self):
        import subprocess,os

        #
        if self.afni_detach:
            os.environ["AFNI_DETACH"] = "YES"
        else:
            os.environ["AFNI_DETACH"] = "NO"


        try:
            if not self.verbose:
                subproc_res = subprocess.run(self.cmd_str, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                
                # log
                if self.log is not None:
                    self.log.write(subproc_res.stdout)
            else:
                # echo out command in formatted way
                self.echo()
                #subproc_res = subprocess.run(self.cmd_str, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                subproc_res = subprocess.run(self.cmd_str, shell=True, universal_newlines=True)
                
                # log
                if self.log is not None:
                    self.log.write("=========\n{}\n=========\n".format(self.cmd_str))
                    self.log.write(subproc_res.stdout)
                    
                #print(subproc_res.stdout)
                
                
            subproc_res.check_returncode()
            
        except subprocess.CalledProcessError as e:
            raise CmdProcError("Error calling external program:\n"+ e.cmd)


    def echo(self):
        """
        print the cmd str to screen in a colorful way reformatted for easy reading
        """
        from termcolor import cprint    #    for colored printing
        cprint( self.cmd_str.replace(" -","\n -").replace("+orig /","+orig\n /").replace(" > ","\n > "), 'green', 'on_grey', attrs=['dark'] )


    def __str__(self):
        """
        override/define print raw cmd str; the 'not-so-colorful' version
        """
        return self.cmd_str


    def append(self,cmd_append=""):
         """
         Wait! i forgot something. Add this to the command.
         """
         self.cmd_str = self.cmd_str + " " + cmd_append


    def set_afni_detach(self,detach):
        """
        By default, if you invoke afni it will NOT detach from the terminal.
        Set this to True if you want to detach.
        """
        self.afni_detach = detach






def prettyOut(msg):
    """
    Standard message output 
    """
    from termcolor import cprint    #    for colored printing

    msg_len = len(msg) + 4
    print("\n")
    cprint("=" * msg_len, 'red', 'on_grey', attrs=['bold'] )
    cprint("= "+ msg + " =", 'red', 'on_grey', attrs=['bold'] )
    cprint("=" * msg_len, 'red', 'on_grey', attrs=['bold'] )
    #print "\n"


def errMsg(msg):

    from termcolor import cprint    #    for colored printing

    msg_len = len(msg) + 4
    print("\n")
    cprint("+" * msg_len, 'yellow', attrs=['bold'] )
    cprint(" ERROR", 'yellow', attrs=['bold'] )
    cprint( " " + msg , 'yellow', attrs=['bold'] )
    cprint("+" * msg_len, 'yellow', attrs=['bold'] )
    #print "\n"



def warningMsg(msg):

    from termcolor import cprint    #    for colored printing

    #msg_len = len(msg) + 4
    print("\n")
    cprint("+++" * 20, 'yellow', attrs=['bold'] )
    cprint(" Warning:", 'yellow', attrs=['bold'] )
    cprint( " " + msg , 'yellow', attrs=['bold'] )
    cprint("+++" * 20, 'yellow', attrs=['bold'] )

    
    

def correct_niftihdr(ref_volume,target_volume,fields):
    """
    Replace field values in target header with those from reference header
    Which fields to replace determined by 'fields' argument
    """
    import subprocess

    proc_res = subprocess.run( args=['nifti_tool','-disp_hdr','-infiles',ref_volume], universal_newlines = True, stdout = subprocess.PIPE)
    proc_res.check_returncode()
    #if proc_res.returncode:
    #    exit(proc_res.returncode)

    # get string return from CompletedProcess object
    ref_hdr = proc_res.stdout
    # build string to write updates to target header
    mod_str    = ""
    for fld_name in fields:
        start_indx     = ref_hdr.find(fld_name)
        end_indx        = ref_hdr.find("\n",start_indx)
        mod_str         = mod_str + " -mod_field " + fld_name + " '" + ' '.join(ref_hdr[start_indx:end_indx].split()[3:]) + "'"

    # write updates to target header
    cmd = ( "nifti_tool -overwrite -infiles " + target_volume + " -mod_hdr " + mod_str )
    return_code = subprocess.call(cmd, shell=True)
    if return_code:
        exit(return_code)


def parse_dicom_path(dicom_path):

    """
    return map of dicom path elements

    """

    path_elements    = dicom_path.split("/")
    seq_elements     = path_elements[-1].split("_")

    directory        = '/'.join(path_elements[:-1])
    sequence         = '_'.join(seq_elements[:-1])
    series           = seq_elements[-1]

    return {'directory': directory, 'sequence': sequence, 'series': series}




def purge_paths(*paths):

    """
    Takes one or more paths inside a list and removes the file objects after expanding any wildcards.

    If path expands into a directory, removes that directory and everything in it.
    If path expands into a file, deletes that file.

    If path has no wildcards, remove that object.

    """

    import glob
    import os
    import shutil


    del_files = []

    #    expand file list with glob
    for p in paths:

        try:
            del_files.extend(glob.glob(p))

        except:
            raise

    #    now purge
    try:
        for fs_obj in del_files:

            if os.path.isdir(fs_obj):
                shutil.rmtree(fs_obj)

            if os.path.isfile(fs_obj):
                os.remove(fs_obj)

    except:
        raise




def make_path(path):
    try:
        os.makedirs(path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def mask_data(i_path, mask_path, o_path):
    import os.path as op
    
    if not op.isfile(o_path):
        try:
            Cmd("fslmaths {} -mas {} {}".format(i_path, mask_path, o_path)).run()
        except:
            raise
            
            
def set_lab_file_permissions(glob_path):
    """
    Quick way to set file permissions to rw for owner and group as you would won't for the lab; removes 'other' permissions
    """
    import glob, stat
    try:
        for f in glob.iglob(glob_path):
            if os.path.isdir(f):
                os.chmod(f, stat.S_IREAD | stat.S_IWRITE | stat.S_IEXEC | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP)
            else:
                os.chmod(f, stat.S_IREAD | stat.S_IWRITE | stat.S_IRGRP | stat.S_IWGRP)
    except:
        raise
            
            
            
            
class AfniImage:
    """
    Class to manipulate AFNI data

    Does not require an AFNI installation.

    SUPPORTED:
        - Supports reading and writing AFNI-like scaled shorts.

    NOT SUPPORTED:
        - different data types for subbriks. Everything in the brik must all be the same data type.
        - subbrik selectors. However, code for this functionality is included but commented out

    UNTESTED:
        - Mac OS or Windows python environment

    """

    # jling Feb 2016.
    #
    # https://afni.nimh.nih.gov/pub/dist/doc/program_help/README.attributes.html


    # this just keeps the import in the class. If AfniImage was moved into its own module you could
    # move this out of this scope
    np = __import__('numpy')


    def __init__(self, data=None, header=None, dtype=None):
        """

        Default constructor initializes only

        Parameters
        ----------

        data
        Optional
        numpy array - becomes copy

        header
        Optional
        dictionary or AFNI attributes - becomes deepcopy

        dtype
        Optional
        numpy dtype
        If not provided, will be determined from AFNI attributes

        Throws TypeError if data is defined but header is not

        """

        if data is not None and header is None:
            raise TypeError("Expected header argument when data argument is defined.")


        self.f_path = None
        self.basename = None
        self.head_file = None
        self.brik_file = None
        self.brik_is_gz = False
        self.data = None
        self.header = None


        if data is not None:

            # data
            if dtype is None:
                self.data = data.copy()
            else:
                self.data = data.astype(dtype, order='F', copy=True)

            # header
            import copy
            self.header = copy.deepcopy(header)




    @classmethod
    def load(cls, filename, dtype=None):
        """

        Return an AfniImage instance

        Parameters
        ----------

        filename
        Required
        Full path to AFNI file.
        Flexible like AFNI binaries. Filename can be 'input+orig' or 'input+orig.HEAD' or 'input+orig.BRIK' or 'input+orig.BRIK.gz', etc.

        dtype
        Optional
        numpy dtype
        If not provided, will be determined from AFNI attributes


        Returns
        -------
        AfniImage

        """

        img = AfniImage()
        img._parse_filename(filename)

        img.header = img._read_hdr_data()

        if dtype is not None:
            img.data = img._read_img_data().astype(dtype)
        else:
            img.data = img._read_img_data()


        return img



    def _parse_filename(self, file_path):
        """

        Provides input filename flexibility.

        Sets this instances path, basename, head and brik filenames

        Throws ValueError if both AFNI files are not found

        Parameters
        ----------
        file_path
            full path and file name

        """

        import os
        self.f_path = os.path.dirname(file_path)
        filename = os.path.basename(file_path)

        test_f = filename.split('.')[-1]
        if (test_f == 'HEAD') or (test_f == 'BRIK'):
            self.basename = '.'.join(filename.split('.')[:-1])
        elif test_f == 'gz':
            self.basename = '.'.join(filename.split('.')[:-2])
        else:
            test_f = filename.split('+')[-1]
            if (test_f == 'orig') or (test_f == 'tlrc'):
                self.basename = filename
            else:
                raise ValueError("Unrecognized file input '{0}'".format(file_path))


        if os.path.isfile( os.path.join(self.f_path,self.basename+".HEAD") ):
            self.head_file = os.path.join(self.f_path,self.basename+".HEAD")
        else:
            raise ValueError("Cannot open HEAD file input '{0}'".format(file_path))

        if os.path.isfile( os.path.join(self.f_path,self.basename+".BRIK") ):
            self.brik_file = os.path.join(self.f_path,self.basename+".BRIK")
        elif os.path.isfile( os.path.join(self.f_path,self.basename+".BRIK.gz") ):
            self.brik_file = os.path.join(self.f_path,self.basename+".BRIK.gz")
            self.brik_is_gz = True
        else:
            raise ValueError("Cannot open BRIK file input '{0}'".format(file_path))




    def _read_hdr_data(self):
        """

        Read the AFNI header file and return dictionary of AFNI attributes

        Returns
        -------
        Dictionary of AFNI attributes

        """

        result = dict()

        with open(self.head_file, "r") as head_file:

            for line in head_file.readlines():

                line_tokens = line.split()

                if not line_tokens:
                    props = { 'type':None, 'count':None, 'name':None, 'value':[] }
                else:
                    if line_tokens[0] in list(props.keys()):
                        (prop, val) = line_tokens[::2]
                        props[prop] = val
                    else:
                        if (props['type'] == 'string-attribute'):
                            # AFNI header strings begin with ' and end with ~
                            props['value'] = line.strip()[1:-1]
                            #props['count'] = len(props['value'])
                        elif (props['type'] == 'float-attribute'):
                            props['value'].extend(list(map(float,line_tokens)))
                        elif (props['type'] == 'integer-attribute'):
                           props['value'].extend(list(map(int,line_tokens)))
                        else:
                            raise ValueError("UnknownProperty")

                        result[props['name']] = props


        return result





    def _read_img_data(self):
        """

        Read and return the raw data from an AFNI BRIK

        Returns
        -------
        ndarray

        """

        # Dictionary translating header field BRICK_TYPES into numpy dtypes.
        # Will throw a KeyError if an unimplemented type is used.
        # JL: changed 0 : "i1" to 0 : "u1" to match AFNI docs
        bricktypes = {
            0 : "u1",
            1 : "i2",
            3 : "f4"
            }

        # Dictionary translating header field BYTEORDER_STRING into numpy dtypes.
        byteorders = {
            "MSB_FIRST" : ">", # Big endian (e.g. PPC Macs)
            "LSB_FIRST" : "<", # Little endian (e.g. x86)
            }


        # check on consistent dtype
        if [a for a in self.header['BRICK_TYPES']['value'] if a != self.header['BRICK_TYPES']['value'][0]]:
            print("ERROR: different types of sub-bricks")
            raise ValueError("DifferentTypes")

        dtype= byteorders[self.header['BYTEORDER_STRING']['value']]
        dtype += bricktypes[self.header['BRICK_TYPES']['value'][0]]

        x,y,z = self.header['DATASET_DIMENSIONS']['value'][0:3]
        nval = self.header['DATASET_RANK']['value'][1]

        try:

            if self.brik_is_gz:
                import gzip
                f=gzip.open(self.brik_file, "rb")
                data=self.np.fromstring(f.read(), dtype=dtype, count=-1)
                f.close()
            else:
                data = self.np.fromfile(file=self.brik_file, dtype=dtype, count=-1)

            # squeeze array shape similar to nibabel data
            if nval == 1:
                data = data.reshape(x, y, z, order='F')
                # check if scaled data
                if any(self.header['BRICK_FLOAT_FACS']['value']):
                    # rescale data and and convert to float
                    data = self.np.array(data*self.header['BRICK_FLOAT_FACS']['value'], dtype='float32', order='F')

            else:
                data = data.reshape(x, y, z, nval, order='F')
                # check if scaled data
                if any(self.header['BRICK_FLOAT_FACS']['value']):
                    # cast array as numpy array so it can be broadcast against multidimensional data
                    float_facs = self.np.array(self.header['BRICK_FLOAT_FACS']['value'])

                    #data = self.np.array(data*self.header['BRICK_FLOAT_FACS']['value'][self.np.newaxis,self.np.newaxis,self.np.newaxis,:], dtype='float32', order='F')
                    # rescale data and and convert to float
                    data = self.np.array(data*float_facs, dtype='float32', order='F')


        except:
            raise

        return data




        # this commented out section provides functionality to load only single
        # or several sub bricks from a volume within the _read_img_data function
        # this is to provide abilities most afni users are accustomed to
        # it might provide improved speed when you only need a single brick from
        # large volume. It could also save memory

        #x = 64
        #y = 64
        #z = 20
        #nval = 90



        #i_file = 'pCASL.fwhm6+orig.BRIK'
        #brik_indx = [9,60]
        #dtype = '<i2'


        ## how to insert gz read?

        #if not brik_indx is None:

            #nbytes = int(dtype[-1])
            #data = self.np.zeros( (x*y*z,len(brik_indx)), dtype=dtype)

            #for indx in range(len(brik_indx)):
                ##f=gzip.open(self.brik_file, "rb")

                #f = open(i_file, "rb")

                #f.seek(x*y*z*nbytes*brik_indx[indx], 0) # size x*y*z*4*1 because 4 bytes/data point
                #data[:,indx] = self.np.fromfile(f, dtype=dtype, count=x*y*z) # count=number of data points
                #f.close()

        #else:
            #data = self.np.fromfile(file=i_file, dtype=dtype, count=-1)



        #if len(brik_indx) == 1:
            #data = data.reshape(x, y, z, order='F')
        #else:
            #data = data.reshape(x, y, z, len(brik_indx), order='F')


        ## now, update header since we only read some bricks





    def _new_data_hdr(self, float_fac=None):
        """

        Return new header AFNI attributes based on current instance to correspond to
        current state of data array. self.data and self.header should be read only on save.
        e.g. shape and dtype

        Parameters
        ----------

        float_fac
        Optional
        If None, data is not scaled and data type determined from self.data
        If not None, data is scaled and header contains necessary info

        Returns
        -------
        new header AFNI attributes based on current instance to correspond to
        current state of data array. self.data and self.header should be read only on save.
        e.g. shape and dtype

        """

        import uuid, copy, os

        # create copy of this instance header
        header = copy.deepcopy(self.header)

        # update header given type and dimension of data
        header['DATASET_DIMENSIONS']['value'][0] = self.data.shape[0]
        header['DATASET_DIMENSIONS']['value'][1] = self.data.shape[1]
        header['DATASET_DIMENSIONS']['value'][2] = self.data.shape[2]

        # i don't want single dimensions in the data shape but AFNI attribute does
        if len(self.data.shape) == 3:
            nval = 1
        else:
            nval = self.data.shape[3]
        header['DATASET_RANK']['value'][1] = nval

        # technically this uuid is too long but afni doesn't mind
        header['IDCODE_STRING']['value'] = "XYZ_{0}".format(uuid.uuid4().hex)

        # this is nice but problematic i'm not sure why i wanted this except that it was cool
        import datetime
        header['HISTORY_NOTE']['value'] = header['HISTORY_NOTE']['value'] + "\\n" + "[{0}@{1}: {2}] {3}".format(os.environ['LOGNAME'],os.uname()[1], datetime.datetime.today().strftime('%c'), 'Save image with AfniImage')


        # BRICK_TYPES and BRICK_FLOAT_FACS
        #
        # create an array of integers that indicates the data type for each subbrik
        # since AfniImage only supports one type for all this creates an array with
        # the length of the time dimension
        if float_fac is not None:

            # set brik data type to int16
            header['BRICK_TYPES']['value'] = self.np.ones((nval,), dtype=self.np.int)
            # include scaling factor in header array
            header['BRICK_FLOAT_FACS']['value'] = self.np.ones((nval,), dtype=self.np.float)*float_fac

        else:

            # set brik data type
            b_out_dtype = self.data.dtype.str[1:3]
            if b_out_dtype == "u1":
                header['BRICK_TYPES']['value'] = self.np.zeros((nval,), dtype=self.np.int)
            elif b_out_dtype == "i2":
                header['BRICK_TYPES']['value'] = self.np.ones((nval,), dtype=self.np.int)
            else:
                header['BRICK_TYPES']['value'] = self.np.ones((nval,), dtype=self.np.int)*3

            # create an array length of time dimension that tells AFNI that data is NOT scaled
            header['BRICK_FLOAT_FACS']['value'] = self.np.zeros((nval,), dtype=self.np.float)



        # BRICK_STATS
        #
        # The main function of this attribute is to provide the display
        # of the dataset numerical ranges on the "Define Overlay" control panel.

        # get min and max for each subbrick
        mn = self.np.nanmin(self.data, axis=(0,1,2))
        mx = self.np.nanmax(self.data, axis=(0,1,2))

        if nval == 1:
            # combine min max for short array
            header['BRICK_STATS']['value'] = [mn, mx]
        else:
            # zip together for a interleaved list of min max for each subbrik
            header['BRICK_STATS']['value'] = [x for mn_mx in zip(mn, mx) for x in mn_mx]


        return header



    def _template_space_filename(self):
        """

        Return proper template space string for filename
        TODO: this needs to be improved
        """
        if self.header['TEMPLATE_SPACE']['value'].lower() != 'orig':
            return 'tlrc'
        else:
            return 'orig'
        
        
    def set_header_value(self, key, value):
        """
        Should properly handle requests for setting or updating header values
        Values can be single values or lists
        
        TODO: why does AFNI use float-attribute for floats and non-floats?
        """
        if key in self.header:
            if self.header[key]['type'] == 'string-attribute':
                self.header[key]['value'] = value
                self.header[key]['count'] = len(value)+1
            elif self.header[key]['type'] == 'float-attribute':
                self.header[key]['value'] = value
                self.header[key]['count'] = len(value)           
        else:
            if isinstance(value, str):
                self.header[key] = {'type': 'string-attribute', 'count': len(value)+1, 'name': key, 'value': value}
            if isinstance(value, list):
                self.header[key] = {'type': 'float-attribute', 'count':len(value), 'name': key, 'value': value }


                
    def save(self, file_path, gz=True, as_short=False):
        """

        Write AfniImage to HEAD/BRIK files

        Parameters
        ----------

        file_path
        Required
        Just the path and what AFNI calls 'prefix'. i.e. Don't include '+orig' or 'HEAD'

        e.g.
        '/tmp/myoutput'
        becomes
        '/tmp/myoutput+orig.HEAD' and '/tmp/myoutput+orig.BRIK.gz' assuming ORIG is the template space

        gz
        Optional
        gzip the output? Defaults to True

        as_short
        Optional
        If True, save float as scaled short AFNI style. Otherwise, float is saved as float.

        """

        import os

        # parse file_path to basename
        f_path = os.path.dirname(file_path)
        basename = os.path.basename(file_path)

        try:

            # user wants to save as AFNI scaled short
            # only do so if self.data is currently a float
            if as_short and self.data.dtype.kind == 'f':

                # get scaling factor by dividing abs(max(data)) by the max of an int16 (32767) to convert to int16
                float_fac = self.np.abs(self.data).max()/32767

                # scale data and cast as int16
                # new array to force correct order
                data = self.np.array(self.data/float_fac, dtype='int16', order='F')

                # update header
                header = self._new_data_hdr(float_fac)

            else:

                # data copy expected

                # AFNI doesn't support float64; if user did some calculation with self.data it could end up
                # as float64.
                # for now, detect and convet
                if self.data.dtype == 'float64':
                    print("Converted output from float64 to float32\n")
                    data = self.data.astype(dtype='float32', order='F', casting='same_kind', copy=True)
                else:
                    data = self.data.copy()

                # sync header before writing to HEAD file
                header = self._new_data_hdr()

        except:
            raise


        try:
            # write brik
            if gz:
                import gzip
                # this method slow on large float datasets. Maybe faster to write unzipped and then gzip separately. Or, switch to native libmri.a
                b_out = gzip.open(os.path.join(f_path, "{0}+{1}.{2}".format(basename, self._template_space_filename(), 'BRIK.gz') ), 'wb')
                b_out.write(data.tobytes(order='F'))
                b_out.close()
            else:
                b_out = open(os.path.join(f_path, "{0}+{1}.{2}".format(basename, self._template_space_filename(), 'BRIK') ), 'wb')
                b_out.write(data.tobytes(order='F'))
                b_out.close()
        except:
            raise




        try:
            # write out the HEAD file
            s_out = open(os.path.join(f_path, "{0}+{1}.{2}".format(basename, self._template_space_filename(), 'HEAD') ), 'w')

            for k, attr in list(header.items()):

                if attr['type'] == 'string-attribute':
                    attr['count'] = len(attr['value'])+1
                    attr['value'] = "'{0}~".format(attr['value'])
                else:
                    attr['count'] = len(attr['value'])
                    attr['value'] = " ".join(map(str,attr['value']))


                # 1 = blank line
                # 2 = type
                s_out.write("\ntype = {0}\n".format(attr['type']))
                # 3 = name
                s_out.write("name = {0}\n".format(attr['name']))
                # 4 = count
                s_out.write("count = {0}\n".format(attr['count']))
                # 5 = the data
                s_out.write("{0}\n".format(attr['value']))

            s_out.close()

        except:
            os.remove(os.path.join(f_path, "{0}+{1}.{2}".format(basename, self._template_space_filename(), 'BRIK.gz') ) )
            raise



    @classmethod
    def concate(cls, file_list, dtype=None):
        """
        A very simple file list to volume converter
        Needs much work
        """
        
        import numpy as np
        
        # how many files in the list
        img_file_count = len(file_list)
        concate_vol_shape = ()
        dim_x, dim_y, dim_z = (0,0,0)
        
        # load each file
        for img_file_indx, img_file_name in enumerate(file_list):
            print(img_file_indx, img_file_name)
            img_file = AfniImage.load(img_file_name)
            
            # at first, get some dimensions and create array
            if not len(concate_vol_shape):
                concate_vol_shape = img_file.data.shape
                dim_x, dim_y, dim_z = img_file.data.shape
                concate_vol = np.zeros( (dim_x, dim_y, dim_z, img_file_count), dtype=img_file.data.dtype)
            elif img_file.data.shape != concate_vol_shape:
                exit("AfniImage:concate:shape shifter!")
            
            # add file data to final volume
            concate_vol[:,:,:,img_file_indx] = img_file.data
            
        # return data volumes
        return concate_vol




class Convert(object):
    """
    A class to contain static methods for conversions
    
    Ungrouped methods are at the top.
    Those methods that make sense to group are put in a subclass to help organize.
    
    """

    ################################################
    # ungrouped methods
    ################################################
    @staticmethod
    def fwhm_to_sigma(fwhm_mm):
        r"""
        Convert FWHM mm to Sigma
        """
        sigma = fwhm_mm/2.35482004503

        return(sigma)
    
    @staticmethod
    def sigma_to_fwhm(sigma):
        """
        Convert Sigma to FWHM 
        """
        fwhm_mm = sigma * 2.35482004503

        return(fwhm_mm)

        

    ################################################
    # grouped methods
    ################################################    
    class Motion(object):

        
        """
        convert rotation degrees to radians
        """
        @staticmethod
        def degreesToRadians(data_degrees):

            from numpy import pi

            # convert degrees to radians ( radians = degrees * pi/180 )
            angle_radians = (data_degrees * (pi/180.0))

            return(angle_radians)
        
        
        
        """
        convert rotation radians to degrees
        """
        @staticmethod
        def radiansToDegrees(data_radians):

            from numpy import pi

            # convert degrees to radians ( radians = degrees * pi/180 )
            angle_degrees = (data_radians * (180.0/pi))

            return(angle_degrees)
        
        
        
        """
        convert rotation degrees to mm
        we will use radius of 50 mm as in JD Powers paper by default
        """
        @staticmethod
        def degreesToMm(data_degrees, radius=50):

            from numpy import pi

            # 1) convert degrees to radians ( radians = degrees * pi/180 )
            angle_radians = (data_degrees * (pi/180.0))

            # 2) using radians and radius, calculate arc length in mm ( arclength = radius * radians(angle))
            arclength = (radius * angle_radians)

            return(arclength)



        """
        convert AFNI conventions to FLIRT and MCFLIRT conventions
        """
        @staticmethod
        def afniToFlirt(in_file_path, out_file_path, skip_lines=0):

            import numpy as np


            # FLIRT and MCFLIRT follow the convention of x_rot y_rot z_rot(all in radians) x_trans y_trans z_trans (all in mm)
            # the parameters are: Rx Ry Rz Tx Ty Tz
            # where R is a rotation (angle in radians) and T is a translation (in mm).

            #Put output in a format similar to 3dvolreg to avoid confusion.
            #i.e.
            #roll pitch yaw dS dL dP
            #
            #                roll  = rotation about the I-S axis }
            #                pitch = rotation about the R-L axis } degrees CCW
            #                yaw   = rotation about the A-P axis }
            #                  dS  = displacement in the Superior direction  }
            #                  dL  = displacement in the Left direction      } mm
            #                  dP  = displacement in the Posterior direction }
            #
            #ASSUMING RPI
            #    z_rot x_rot y_rot z_trans x_trans y_trans
            #    1     2     3     4       5       6
            #
            #
            #R = Y X Z where Z is rotation about the I-S axis, X is rotation about the R-L axis, and Y is rotation about the A-P axis


            # ORDER: import 3dVolreg output columns in an order that matches FLIRT order
            mc_params = np.genfromtxt(in_file_path, usecols=(2,3,1,5,6,4), dtype=np.float64, skip_header=skip_lines)

            # UNIT: convert rot on 50mm sphere
            mc_params[:,0:3] = Convert.Motion.degreesToRadians(mc_params[:,0:3])

            # ORIENT: fix sign on columns 2(1) and 5(4) to match MCFLIRT
            # i'm assuming -orient RPI
            # this would make it exactly match output of something like MCFLIRT but here i didn't want to do that
            #mc_params[:,(1,4)] *= -1

            np.savetxt(out_file_path, mc_params, fmt="%.5e", delimiter="  ")



class Gpu(object):
        
        r""" 
        Returns index of unused or underused gpu id
        Used to specify GPU in environment variables
        """
        @staticmethod
        def get_available(util_threshold=25.0):

            import nvsmi
            
            # return unused gpu                
            for gpu in nvsmi.get_available_gpus():
                return(gpu.id)

            # if they are utilized, get by threshold
            for gpu in nvsmi.get_gpus():
                if gpu.gpu_util < util_threshold:
                    return gpu.id
                
            return None
        
        
        r""" 
        Returns number of available gpus
        """
        @staticmethod
        def get_gpu_count():
            import nvsmi
            try:
                gpu_count = len(list(nvsmi.get_gpus()))
            except:
                gpu_count = 0
            finally:
                return gpu_count
        
        
        
        @staticmethod
        def get_gpu_names():
            import nvsmi
            for gpu in nvsmi.get_gpus():
                return gpu.name

        


        
