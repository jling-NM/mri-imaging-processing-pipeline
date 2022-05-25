#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""@package anat

    Created on Wed Mar  3 12:08:48 2021
    
    @author: jling

    PRISMA
    
    Processing module for anat data.
    
    BIDS-dependent
 
"""

import os
import errno
from mayerlab.preprocessing import utils


def cortical_thickness(subj, study, template_name='tpl.v1_v2'):
    """
    Generate ANTs input for mindboggle
    
    antsCorticalThickness.sh performs T1 anatomical brain processing where the following steps are currently applied:
      1. Brain extraction
      2. Brain n-tissue segmentation
      3. Cortical thickness
      4. (Optional) registration to a template
              
    """
    
    template_path = '/export/research/analysis/human/amayer/shared/apps/mindboggle/OASIS-30_Atropos_template'


    try:
        
        visit_num = str(int(subj.visit[-1]))
        
        long_dir = os.path.join(study.freesurfer_path, 'longitudinal', subj.ursi.full, 'v' + visit_num + '.long.' + template_name)
        ants_dir = os.path.join(long_dir, 'ants')

        if not os.path.isfile(os.path.join(ants_dir, 'ACTStage6Complete.txt' )):
            
            ants_input = os.path.join(long_dir, 'mri', 'T1.mgz')

            if not os.path.isfile(ants_input):
                raise ValueError('mindboggle_ants_cortical_thickness() expects processed freesurfer output at ' + ants_input)
                
            # make ants path
            try:
                os.makedirs(ants_dir, 0o770)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            
            os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = '6'

            # run ants
            utils.Cmd('antsCorticalThickness.sh -d 3 -a ' + ants_input +
                ' -o ' + ants_dir + '/' +  
                ' -e ' + os.path.join(template_path, "T_template0.nii.gz") + 
                ' -t ' + os.path.join(template_path, "T_template0_BrainCerebellum.nii.gz") + 
                ' -m ' + os.path.join(template_path, "T_template0_BrainCerebellumProbabilityMask.nii.gz") + 
                ' -f ' + os.path.join(template_path, "T_template0_BrainCerebellumExtractionMask.nii.gz") + 
                ' -p ' + os.path.join(template_path, "Priors2/priors%d.nii.gz") + 
                ' -u 0').run()

    except:
        raise
        
    
    
    
    
def dkt_atlas(subj, study, template_name='tpl.v1_v2'):
    """        
    # add DKT31atlas101 template to longitudinal processed data for mindboggle
    # lh.DKT31atlas101subjects.gcs # OSF mindboggle 101
    # lh.DKTatlas100.gcs # fs 6.0
      
    """
    
    from mayerlab import utils
    import glob, shutil
    
    gcs_i_file_sufx = "DKT31atlas101subjects"
    gcs_o_file_sufx = "DKTatlas101"


    try:
        
        os.environ["OMP_NUM_THREADS"] = '6'
                    
        visit_num = str(int(subj.visit[-1]))
        
        sdir     = os.path.join(subj.ursi.full, 'v' + visit_num + '.long.' + template_name)
        long_dir = os.path.join(study.freesurfer_path, 'longitudinal', sdir)
        
        for hemi in ['lh', 'rh']:

            # create DKTatlas101 annot file
            if not os.path.isfile(os.path.join(long_dir, 'label', 'annot', hemi+'.aparc.'+gcs_o_file_sufx+'.annot')):

                utils.Cmd("mris_ca_label -seed 1234 -sdir {} -l {} -aseg {} {} {} {} {} {}".format(
                    os.path.join(study.freesurfer_path, 'longitudinal'),
                    os.path.join(long_dir, 'label', hemi + '.cortex.label'),
                    os.path.join(long_dir, 'mri', 'aseg.mgz'),
                    sdir,
                    hemi,
                    os.path.join(long_dir, 'surf', hemi + '.sphere.reg'),
                    os.environ['FREESURFER_HOME'] + '/average/' + hemi + '.' + gcs_i_file_sufx + '.gcs',
                    os.path.join(long_dir, 'label', hemi + '.aparc.' + gcs_o_file_sufx + '.annot')
                    )).run()
                
                    
            
                # create DKTatlas101 ctab
                i_file = os.path.join(long_dir, 'label', hemi+'.aparc.'+gcs_o_file_sufx+'.annot')
                o_file_ctab = os.path.join(long_dir, 'label', 'aparc.annot.'+gcs_o_file_sufx+'.ctab')
                o_file_stats = os.path.join(long_dir, 'stats', hemi+'.aparc.'+gcs_o_file_sufx+'.stats')
                if not os.path.isfile(o_file_stats):
                    utils.Cmd("mris_anatomical_stats -sdir {} -mgz -cortex {} -f {} -b -a {} -c {} {} {} white".format(
                        os.path.join(study.freesurfer_path, 'longitudinal'),
                        os.path.join(long_dir, 'label', hemi + '.cortex.label'),
                        o_file_stats,
                        i_file,
                        o_file_ctab,
                        sdir,
                        hemi
                        )).run()



        # create DKTatlas101 parcellation volume
        if not os.path.isfile(os.path.join(long_dir, 'mri', 'aparc', 'aparc.'+gcs_o_file_sufx+'+aseg.mgz')):
            
            os.environ['SUBJECTS_DIR'] = os.path.join(study.freesurfer_path, 'longitudinal')
            
            utils.Cmd("mri_aparc2aseg --s {} --volmask --annot {}".format(
                sdir,
                'aparc.'+gcs_o_file_sufx
                )).run()
            
        
        # create DKTatlas101 wmparc volume
        o_file = "{}/mri/wmparc/wmparc.{}.mgz".format(long_dir, gcs_o_file_sufx)
        
        if not os.path.isfile(o_file):
            os.environ['SUBJECTS_DIR'] = os.path.join(study.freesurfer_path, 'longitudinal')
            
            utils.Cmd("mri_aparc2aseg --s {} --labelwm --hypo-as-wm --rip-unknown --volmask --annot aparc.DKTatlas101 --o {} --ctxseg {}".format(
            sdir,
            "{}/mri/wmparc.{}.mgz".format(long_dir, gcs_o_file_sufx),
            'aparc.'+gcs_o_file_sufx+'+aseg.mgz'
            )).run()
            

        
        # since mindboggle uses default freesurfer file inputs we need to move stuff and create pointers so mindboggle uses the correct DKT101 files
        
        dest_dir = os.path.join(long_dir, 'mri', 'aparc')
        if not os.path.isdir(dest_dir):
            try:
                os.makedirs(dest_dir, 0o770)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

            for f in glob.iglob(os.path.join(long_dir, 'mri', 'aparc*aseg.mgz')):
                shutil.move(f, os.path.join(long_dir, 'mri', 'aparc'))
        
            os.symlink( os.path.join(long_dir, 'mri', 'aparc', 'aparc.'+gcs_o_file_sufx+'+aseg.mgz'), os.path.join(long_dir, 'mri', 'aparc+aseg.mgz') )
        
        
        
        
        
        dest_dir = os.path.join(long_dir, 'mri', 'wmparc')
        if not os.path.isdir(dest_dir):
            try:
                os.makedirs(dest_dir, 0o770)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
                    
            for f in glob.iglob(os.path.join(long_dir, 'mri', 'wmparc*')):
                shutil.move(f, os.path.join(long_dir, 'mri', 'wmparc'))
        
            os.symlink( os.path.join(long_dir, 'mri', 'wmparc', 'wmparc.'+gcs_o_file_sufx+'.mgz'), os.path.join(long_dir, 'mri', 'wmparc.mgz') )
        
        
        
        
        dest_dir = os.path.join(long_dir, 'label', 'annot')
        if not os.path.isdir(dest_dir):
            try:
                os.makedirs(dest_dir, 0o770)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

            for f in glob.iglob(os.path.join(long_dir, 'label', '?h.aparc.*')):
                shutil.move(f, os.path.join(long_dir, 'label', 'annot'))
        
            for hemi in ['lh', 'rh']:
                os.symlink( os.path.join(long_dir, 'label', 'annot', hemi+'.aparc.'+gcs_o_file_sufx+'.annot'), os.path.join(long_dir, 'label', hemi+'.aparc.annot') )

        
        
        dest_dir = os.path.join(long_dir, 'label', 'ctab')
        if not os.path.isdir(dest_dir):
            try:
                os.makedirs(dest_dir, 0o770)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
                    
            for f in glob.iglob(os.path.join(long_dir, 'label', 'aparc*.ctab')):
                shutil.move(f, os.path.join(long_dir, 'label', 'ctab'))
                
            os.symlink( os.path.join(long_dir, 'label', 'ctab', 'aparc.annot.'+gcs_o_file_sufx+'.ctab'), os.path.join(long_dir, 'label', 'aparc.annot.ctab') )


    except:
        raise
    
    
    
        
def estimates(subj, study, template_name='tpl.v1_v2'):
    
    from mayerlab import utils
    
    try:
        
        os.environ["OMP_NUM_THREADS"] = '6'
                    
        visit_num = str(int(subj.visit[-1]))
        
        sdir     = os.path.join(subj.ursi.full, 'v' + visit_num + '.long.' + template_name)
        long_dir = os.path.join(study.freesurfer_path, 'longitudinal', sdir)            
        ants_dir = os.path.join(long_dir, 'ants')
        boggled_dir = os.path.join(long_dir, 'boggled')

        if not os.path.isfile( os.path.join(boggled_dir, 'tables', 'left_cortical_surface', 'fundus_shapes.csv') ):
            
            try:
                os.makedirs(boggled_dir, 0o770)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
                    
            os.environ["vtk_cpp_tools"] = '/export/research/analysis/human/amayer/shared/apps/mindboggle/vtk_cpp_tools'
            os.environ["LD_LIBRARY_PATH"] = '/export/research/analysis/human/amayer/shared/apps/VTK/lib:/export/research/analysis/human/amayer/shared/apps/anaconda/lib'
            
            # longitudinal data missing contents of orig folder that mindboggle expects
            # it is not entirely clear what the purpose of this is and i get no response from postings
            # Freesurfer docs say: "After the registrations and the norm_template.mgz volume are created, the orig.mgz images from all TPs are mapped to the template location and averaged (again the default is the median image) to produce the <templateid>/mri/orig/001.mgz image"
            # but i don't see that file
            # the pipeline does provide an orig.mgz and i will link that for now. It is already in the same coordinate system as the T1.mgz
            target_link_path = os.path.join(long_dir, 'mri', 'orig', '001.mgz')
            if not os.path.isfile(target_link_path):
                os.symlink(os.path.join(long_dir, 'mri', 'orig.mgz'), target_link_path)
            
            # mindboggle on both freesurfer and ANTs
            # note: using "--roygbiv" crashes out
            #mindboggle ${long_src} --out ${boggled_dir} --ants ${ants_dir}/BrainSegmentation.nii.gz --cpus 6 --fundi --working /tmp/mindboggle | tee -a ${boggled_dir}/mindboggle.run.log
            utils.Cmd('mindboggle '+long_dir+' --out '+boggled_dir+' --ants '+ants_dir+'/BrainSegmentation.nii.gz --cpus 6 --fundi --working /tmp/mindboggle/'+subj.ursi.full+' | tee -a '+boggled_dir+'/mindboggle.run.log').run()
            
    except:
        raise
                
