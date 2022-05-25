#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 11:37:16 2021

@author: jling
"""



def standard(subj):
    """
    New lab standard. Currently used for MCTBI
    """
    
    from mayerlab.preprocessing import utils
    from mayerlab.preprocessing.anat import afni
    from mayerlab.preprocessing.anat import ants
    from mayerlab.preprocessing.anat import freesurfer
    from mayerlab.preprocessing.anat import spm
    from mayerlab.preprocessing.anat import qc
    
    utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT Processing")
    
    # crop FOV
    afni.crop_t1(subj)
    
    # brain mask
    ants.brain_extraction(subj)

    # T1 to MNI
    ants.warp_t1_to_MNI(subj, template_key='mni_icbm152_nlin_asym_2009c')
    qc.snapshot_t1_to_MNI(subj)
    
    # spm priors tissue segmentation
    spm.generate_tissue_classes(subj)
    
    qc.run_mriqc(subj)
    
    freesurfer.recon_all(subj, use_T2=True)

    qc.aggregate_qc_pdf(subj)
    
    utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT Processsing Complete")
    
    
    
    
    
    
def laptop(subj):
    """
    Pipeline to closely replicate processsing used for TRIO data on early Laptop data
    For the exact Trio pipeline see /mayerlab/preprocessing/trio
    """
    from mayerlab.preprocessing import utils
    from mayerlab.preprocessing.anat import afni
    from mayerlab.preprocessing.anat import freesurfer
    from mayerlab.preprocessing.anat import spm
    from mayerlab.preprocessing.anat import qc



    utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT Processsing")
    
    
    # these functions duplicate laptop_20753 spatial normalization for laptop_20019
    afni.skullstrip(subj)
    afni.warp_t1_to_TLRC(subj)
    qc.snapshot_t1_to_TLRC(subj)
    
    
    # run freesurfer
    freesurfer.recon_all(subj, use_T2=True)
    # T1 segmentation
    spm.generate_tissue_classes(subj)
    # mriqc 
    qc.run_mriqc(subj)     
    # put together qc images
    qc.aggregate_qc_pdf(subj, space='TLRC')
    
    utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT Processsing Complete")

