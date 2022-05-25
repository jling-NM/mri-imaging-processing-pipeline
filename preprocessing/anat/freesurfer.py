#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""@package anat

    Created on Wed Mar  3 12:08:48 2021
    
    @author: jling

    PRISMA
    
    Processing module for anat data.
    
    BIDS-dependent
 
"""

from mayerlab.preprocessing import utils
from mayerlab.preprocessing.anat import afni
import os




def recon_all(subj, flag_mprage=True, use_T2=False):
    
    """
    Runs FreeSurfer recon-all on T1.
    
    :param subj: utils.Subj
    :param study: utils.Study    
    :param flag_mprage: Should recon-all run with "-mprage" flag?
    
    
    
    Using T2 or FLAIR data to improve pial surfaces
    Pial surfaces can be improved using the different contrast in T2 or FLAIR images. The original pial surfaces without T2/FLAIR data are saved as ?h.woT2.pial or ?h.woFLAIR.pial, and new ?h.pial surfaces are created. One example where this is useful is when there is dura in the brainmask.mgz that isn't removed by skull stripping. The flags for these commands are:
    -T2 <input T2 volume> or -FLAIR <input FLAIR volume>
    	
    
    specify the path to the T2 or FLAIR image to use (can be either a DICOM, MGH or NIFTI file)
    -T2pial or -FLAIRpial
    	
    creates new pial surfaces using T2 or FLAIR images
    
    An example of running a subject through Freesurfer with a T2 image is:
    recon-all -subject subjectname -i /path/to/input_volume -T2 /path/to/T2_volume -T2pial -all
    
    T2 or FLAIR images can also be used with Freesurfer subjects that have already been processed without them. Note that autorecon3 should also be re-ran to compute statistics based on the new surfaces. For example:
    recon-all -subject subjectname -T2 /path/to/T2_volume -T2pial -autorecon3

    
    
    """
    
    try:
    
        # we want to crop the T1 before using for anything; probably already
        afni.crop_t1(subj)
        
        anat_path = subj.derivatives_path(bidsType='anat')      
        # this needs to be here before first subject is run
        utils.make_path(subj.study.freesurferPath)
    
        # freesurfer looks in environment for this
        os.environ['SUBJECTS_DIR'] = subj.study.freesurferPath
        
        # check for fs output first
        if os.path.exists(os.path.join(subj.study.freesurferPath, subj.ursi.full)):

            # check for partial completion
            if not os.path.isfile(os.path.join(subj.study.freesurferPath,subj.ursi.full,"mri","brainmask.mgz")) and not os.path.isfile(os.path.join(subj.study.freesurferPath,subj.ursi.full,"mri","T1.mgz")):
                exit("Half-baked Freesurfer output?")
            
            else:
                utils.prettyOut("Freesurfer recon_all output already exists")
                
        else:

            t1_in_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
            if not os.path.isfile(t1_in_path):
                raise utils.CmdProcError("T1w input not found")
                

            # run recon_all
            utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT : FreeSurfer" )

            cmd = utils.Cmd("recon-all -all -parallel -threads 9 -3T -i "+ t1_in_path +" -subject "+subj.ursi.full)
            
            if flag_mprage:
                cmd.append("-mprage")
                
            if use_T2:
                t2_in_path = subj.bids_file_path(bidsType='anat', bidsLabel='T2w')
 
                if t2_in_path != None:
                    cmd.append("-T2 "+ t2_in_path +" -T2pial")
                else:
                    utils.warningMsg("Requested T2w input not found. \n FreeSurfer proceeding with T1w.")
                
            cmd.run()

                        
    except:
        
        raise
    
    
    
    

def gen_qc_json(subj, study):
    """
    Create a single-subject summary of etiV normalized volumes, and average thickness in json format
    """
    
    import pandas as pd
    
    try:
        
        subj_path = os.path.join(study.freesurfer_path,  subj.visit, subj.ursi.full, 'stats')
        
        # summarize all subcortical volumes
        input_path = os.path.join(subj_path, 'aseg.stats')
        output_path = os.path.join(subj_path, 'aseg.json')
        
        if not os.path.isfile(output_path):
            
            with open(input_path, 'r') as f:
                for line in f.readlines():
                    if 'EstimatedTotalIntraCranialVol' in line:
                        eTIV = float(line.split(',')[3].strip())
                        break
        
            df = pd.read_table(input_path, 
                               comment='#', 
                               sep='\s+', 
                               names=['Index', 'SegId', 'NVoxels', 'Volume_mm3', 'StructName', 'normMean', 'normStdDev', 'normMin', 'normMax', 'normRange'],
                               index_col='StructName'
                              )
            df['Volume_eTIV_norm'] = df['Volume_mm3']/eTIV
        
            df.to_json(output_path, orient="index", indent=4)
            print("   +++ Wrote file:", output_path)
        
        
        
        
        # summarize all cortical GM
        for hemi in ['lh','rh']:
            
            input_path = os.path.join(subj_path, hemi+'.aparc.stats')
            output_path = os.path.join(subj_path, hemi+'.aparc.json')
        
            if not os.path.isfile(output_path):
                with open(input_path, 'r') as f:
                    for line in f.readlines():
                        if 'EstimatedTotalIntraCranialVol' in line:
                            eTIV = float(line.split(',')[3].strip())
                            break
        
                df = pd.read_table(input_path, 
                                   comment='#', 
                                   sep='\s+', 
                                   names=['StructName','NumVert','SurfArea','GrayVol','ThickAvg','ThickStd','MeanCurv','GausCurv','FoldInd','CurvInd'],
                                   index_col='StructName'
                                  )
                df['GrayVol_eTIV_norm'] = df['GrayVol']/eTIV
        
                df.to_json(output_path, orient="index", indent=4)
                print("   +++ Wrote file:", output_path)

    except:
        raise
        
        
        
def get_qc_json(subj, study):
    """
    Returns joined stats in json format
    """
    
    import json
    
    try:
        
        subj_path = os.path.join(study.freesurfer_path,  subj.visit, subj.ursi.full, 'stats')
        
        # start master dictionary here
        all_qc = {'URSI':subj.ursi.full, 'VISIT':subj.visit, 'STATS': {'aseg':None, 'lh.aparc':None, 'rh.aparc':None} }
        
        # get each item and add to master
        for (region, _) in all_qc['STATS'].items():
            
            input_path = os.path.join(subj_path, region+'.json')
            if not os.path.isfile(input_path):
                gen_qc_json(subj, study)
                
            with open(input_path,'r') as f:
                all_qc['STATS'][region] = json.load(f)

            
        # return joined master
        #return json.dumps(all_qc, indent = 4)
        return all_qc


    except:
        raise
        
    
    

    
    
def xhemi(subj,study):
    """
    Surface-based Interhemispheric Registration
    https://surfer.nmr.mgh.harvard.edu/fswiki/Xhemi
    
    Greve, Douglas N., Lise Van der Haegen, Qing Cai, Steven Stufflebeam, Mert R. Sabuncu, Bruce Fischl, and Marc Bysbaert. "A surface-based analysis of language lateralization and cortical asymmetry." (2013). Journal of Cognitive Neuroscience 25.9: 1477-1492.
    
    """
 
        
    # freesurfer needs link for fsaverage_sym in SUBJECTS_DIR
    fs_sym_path = os.path.join(study.freesurfer_path, subj.visit, 'fsaverage_sym')        
    if not os.path.islink(fs_sym_path):
        os.symlink(os.path.join(utils.Env.freesurfer_path, 'subjects', 'fsaverage_sym'), fs_sym_path)
        
    #    let's not bombard hosts
    os.environ['OMP_NUM_THREADS'] = '6'
        
    # Reg to atlas
    if not os.path.isfile(os.path.join(subj.get_path(study, 'fs'), 'surf', 'lh.fsaverage_sym.sphere.reg')):
        # Creates lh.fsaverage_sym.sphere.reg in $subject and $subject/xhemi
        utils.prettyOut("{} : {} : Surface-based Interhemispheric Registration : Stage 1".format(subj.visit, subj.ursi.full))
        utils.Cmd("surfreg --s "+ subj.ursi.full + " --t fsaverage_sym --lh").run()
        
    if not os.path.isfile(os.path.join(subj.get_path(study, 'fs'), 'xhemi', 'lrrev.register.dat')):
        utils.prettyOut("Surface-based Interhemispheric Registration : {} : {} : Stage 2".format(subj.visit, subj.ursi.full))
        utils.Cmd("surfreg --s "+ subj.ursi.full + " --t fsaverage_sym --lh --xhemi").run()
