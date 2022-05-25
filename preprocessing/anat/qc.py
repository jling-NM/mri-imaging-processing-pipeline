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
from mayerlab.preprocessing.anat import ants
import os



def run_mriqc(subj):
    """
    
    """

    if ( (not os.path.isfile(os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'anat', subj.bids_file_join('T1w.json'))))
        or (not os.path.isfile(os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'anat', subj.bids_file_join('T2w.json'))))
    ):

        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT : MRIQC")

        # must have this env variable set to have access to templates; does not work to include as --env 
        os.environ['SINGULARITYENV_TEMPLATEFLOW_HOME'] = "/home/bidsapp/.cache/templateflow"
        utils.Cmd("singularity run --cleanenv --containall --writable-tmpfs --bind /export:/export /export/research/analysis/human/amayer/shared/apps/containers/runtime/mayerlab_mriqc_0.16.1.sif "
                  + " --participant_label " + subj.ursi.short 
                  + " --session-id " + str(subj.visit)
                  + " --modalities T1w T2w "
                  + " --work-dir " + subj.study.mriqcPath + "/workspace/" + subj.bids_file_join()
                  + " --no-sub "
                  + " --verbose "
                  + " --nprocs 6 "
                  + " --verbose-reports "
                  + " --float32 "
                  + " --ants-nthreads 6 "
                  + " --omp-nthreads 6 "
                  + subj.study.bidsPath + " "
                  + subj.study.mriqcPath + " participant "
                  ).run()
        
        
        # have to clean up after poldrack
        utils.purge_paths(subj.study.mriqcPath + "/workspace/" + subj.bids_file_join())
        
        try:
            # fix permissions on container task output for group access
            utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id) )
            utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit) )
            utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_file_join('T1w.html')) )
            utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_file_join('T2w.html')) )
            utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'anat') )
            utils.set_lab_file_permissions( os.path.join(subj.study.mriqcPath, subj.bids_id, subj.bids_visit, 'anat', '*.json') )
        except:
            pass        





def snapshot_tissue_classes(subj):
    
    
    """
    
    """    
    
    anat_path = subj.derivatives_path(bidsType='anat')
    ulay_file = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))

    utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT : QC Segmentation Snapshots" )

    tissue_dict = dict(gm='1',wm='2',csf='3')
    for key in tissue_dict:
        mask_file = os.path.join(anat_path, subj.bids_file_join('spm-' + key + '_mask.nii.gz'))
        
        utils.Cmd("@chauffeur_afni \
                -ulay "+ulay_file+" \
                -olay "+mask_file+" \
                -prefix "+ os.path.join(anat_path, 'qc', subj.bids_file_join('spm-' + key + '_mask')) +" \
                -set_xhairs OFF \
                -label_string " + subj.bids_file_join('spm-' + key + '_mask') + " \
                -montx 4 \
                -monty 1 \
                -label_string "+subj.bids_file_join('spm-' + key + '_mask') +" \
                -delta_slices 16 16 16 \
                -pbar_posonly \
                -cbar_ncolors 1").run()
    




def snapshot_t1_to_TLRC(subj):
    
    
    """
    
    """    
    
    
    
    anat_path = subj.derivatives_path(bidsType='anat')
    # make sure qc subdirectory exists
    utils.make_path(os.path.join(anat_path, 'qc'))

    test_file = os.path.join(anat_path, 'qc', subj.bids_file_join('T1w-TLRC.jpg'))
    if not os.path.isfile(test_file):
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT : QC T1w to TLRC Snapshots" )
        
        ulay_file = utils.Env.afni_path + "/TT_N27+tlrc"
        olay_file = os.path.join(anat_path, subj.bids_file_join('T1w_SKSP.TLRC.nii.gz'))

        utils.Cmd("@snapshot_volreg {} {} {}".format(ulay_file, olay_file, os.path.join(anat_path, 'qc', subj.bids_file_join('T1w-TLRC')))).run()






def snapshot_t1_to_MNI(subj):
    """
    """    
    
    anat_path = subj.derivatives_path(bidsType='anat')
    # make sure qc subdirectory exists
    utils.make_path(os.path.join(anat_path, 'qc'))

    test_file = os.path.join(anat_path, 'qc', subj.bids_file_join('T1w-MNI.jpg'))
    if not os.path.isfile(test_file):
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT : QC T1w to MNI Snapshots" )
        
        
        i_path = os.path.join(anat_path, subj.bids_file_join('T1w.nii.gz'))
        o_path = os.path.join(anat_path, subj.bids_file_join('VERIFY.T1w-to-MNI.nii.gz'))
        
        ants.apply_spatial_normalization(subj, template_key='mni_icbm152_nlin_asym_2009c', 
                                         i_path=i_path, 
                                         o_path=o_path)
        
        
        ulay_file = utils.Env.template_lib['mni_icbm152_nlin_asym_2009c']['head']
        olay_file = o_path
        
        utils.Cmd("@snapshot_volreg {} {} {}".format(ulay_file, olay_file, os.path.join(anat_path, 'qc', subj.bids_file_join('T1w-MNI')))).run()

        utils.purge_paths(o_path)
        
        
        
        
        
        
def aggregate_qc_pdf(subj, space='MNI'):   
    """
    
    """    
    import pylatex
    from pylatex.utils import NoEscape        
    
    
    anat_path = subj.derivatives_path(bidsType='anat')
    
    o_file = os.path.join(anat_path, 'qc', 'qc-preprocessing')
    if not os.path.isfile(o_file + '.pdf'):
        
        utils.prettyOut(subj.ursi.full + " : Visit " + str(subj.visit) + " : ANAT : QC PDF" )
        
        
        geometry_options = {"tmargin": "1cm", "lmargin": "3cm"}
        doc = pylatex.Document(page_numbers=False, geometry_options=geometry_options)

        doc.preamble.append(pylatex.Command('title', subj.study.label + ' ANAT Preprocessing QC'))
        doc.preamble.append(pylatex.Command('author', subj.ursi.full + " : Visit " + str(subj.visit)))
        doc.preamble.append(pylatex.Command('date', NoEscape(r'\today')))
        doc.append(NoEscape(r'\maketitle'))
        doc.append(NoEscape(r'\newpage'))
        
        
        with doc.create(pylatex.Center()):
            with doc.create(pylatex.Section('Spatial Normalization', numbering=False)):
                doc.append('Colored edge overlay should align with tissue boundaries')
                with doc.create(pylatex.Figure(position='h')) as pic1:
                    pic1.add_image( os.path.join(anat_path, 'qc', subj.bids_file_join('T1w-'+space+'.jpg')), width=NoEscape('1.0\\linewidth'), placement=None )
                    pic1.add_caption('Subject T1w registration over '+space+' template')
                
        doc.append(NoEscape(r'\newpage'))
        
        with doc.create(pylatex.Center()):
            with doc.create(pylatex.Section('Tissue Masks', numbering=False)):
                with doc.create(pylatex.Figure(position='h')) as pic1:
                    pic1.add_image(os.path.join(anat_path, 'qc', subj.bids_file_join('spm-gm_mask.axi.png')), width=NoEscape('0.8\\linewidth'), placement=NoEscape('\\centering') )
                    pic1.add_caption('GM Mask - Axial View')
                with doc.create(pylatex.Figure(position='h')) as pic2:
                    pic2.add_image(os.path.join(anat_path, 'qc', subj.bids_file_join('spm-wm_mask.axi.png')), width=NoEscape('0.8\\linewidth'), placement=NoEscape('\\centering') )
                    pic2.add_caption('WM Mask - Axial View')
                with doc.create(pylatex.Figure(position='h')) as pic3:
                    pic3.add_image(os.path.join(anat_path, 'qc', subj.bids_file_join('spm-csf_mask.axi.png')), width=NoEscape('0.8\\linewidth'), placement=NoEscape('\\centering') )
                    pic3.add_caption('CSF Mask - Axial View')  
                
        doc.append(NoEscape(r'\newpage'))
        
        with doc.create(pylatex.Center()):
            with doc.create(pylatex.Section('Tissue Masks', numbering=False)):
                with doc.create(pylatex.Figure(position='h')) as pic1:
                    pic1.add_image(os.path.join(anat_path, 'qc', subj.bids_file_join('spm-gm_mask.sag.png')), width=NoEscape('0.8\\linewidth'), placement=NoEscape('\\centering') )
                    pic1.add_caption('GM Mask - Sagittal View')
                with doc.create(pylatex.Figure(position='h')) as pic2:
                    pic2.add_image(os.path.join(anat_path, 'qc', subj.bids_file_join('spm-wm_mask.sag.png')), width=NoEscape('0.8\\linewidth'), placement=NoEscape('\\centering') )
                    pic2.add_caption('WM Mask - Sagittal View')
                with doc.create(pylatex.Figure(position='h')) as pic3:
                    pic3.add_image(os.path.join(anat_path, 'qc', subj.bids_file_join('spm-csf_mask.sag.png')), width=NoEscape('0.8\\linewidth'), placement=NoEscape('\\centering') )
                    pic3.add_caption('CSF Mask - Sagittal View')  
                    
                    
        doc.generate_pdf(o_file, clean_tex=True, compiler='pdflatex')
        
        
        