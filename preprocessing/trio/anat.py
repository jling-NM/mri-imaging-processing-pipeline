"""@package anat
Preprocessing module for T1
 
"""

from mayerlab import utils    #
from sys import exit        #    exit()
import os                #
import subprocess            #    subprocess.check_call()
import errno            #    for path check and create
import glob



      
    
    
def mindboggle_ants_cortical_thickness(subj, study):
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
        
        ants_dir = os.path.join(study.path, subj.ursi.full, subj.visit, 'anat', 'ants')

        if not os.path.isfile(os.path.join(ants_dir, 'ACTStage6Complete.txt' )):
            
            ants_input = os.path.join(study.freesurfer_path, subj.visit, subj.ursi.full, 'mri', 'T1.mgz')

            if not os.path.isfile(ants_input):
                raise ValueError('mindboggle_ants_cortical_thickness() expects processed freesurfer output at ' + ants_input)
                
            # make ants path
            try:
                os.makedirs(ants_dir, 0o770)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            
            # run ants
            os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = '6'
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
        
        
        
        

def gen_freesurfer_qc_json(subj, study):
    """
    Create a single-subject summary of etiV normalized volumes, and average thickness in json format
    """
    
    import os
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
        
        
        
def get_freesurfer_qc_json(subj, study):
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
                gen_freesurfer_qc_json(subj, study)
                
            with open(input_path,'r') as f:
                all_qc['STATS'][region] = json.load(f)

            
        # return joined master
        #return json.dumps(all_qc, indent = 4)
        return all_qc


    except:
        raise
        
        
        
def to3d(subj,study):
    
    """
    Volumize dicom input, crops image
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    out_path = study.path + "/" + subj.ursi.full + subj.visit_path + "/anat/"
    out_file = study.label + "." + subj.ursi.short + ".anat"
    
    #    try to make output path and fail for reasons other than it already exists
    try:
        os.makedirs(out_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

        
    if not os.path.isfile(out_path + out_file + "+orig.HEAD"):

        #    need dicom path from _run_map_anat_
        dcm_path = subj.get_dcm_path(study,"anat","MPRAGE1")
        if not dcm_path:
            raise utils.CmdProcError("Dicom path not found")
            
        utils.prettyOut("Create anat volume...")
        
        utils.Cmd( "to3d -anat -session " + out_path + " -prefix anat.tmp " + dcm_path + "/*.dcm" ).run()
        if os.path.isfile(os.path.join(out_path,"anat.tmp+orig.HEAD")):
                    utils.Cmd("3dautobox -input " + os.path.join(out_path,'anat.tmp+orig.HEAD') +" -prefix " + os.path.join(out_path, out_file) ).run()
                    utils.purge_paths(os.path.join(out_path,"anat.tmp+orig.*"))
            
            
            

def assemble_volume(subj,study):
    
    
    """
    Robust volumizing of dicom input.
    Good foreign data, non Siemens data as it more carefully checks image count using the header.        
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    out_path = study.path + "/" + subj.ursi.full + subj.visit_path + "/anat/"
    out_file = study.label + "." + subj.ursi.short + ".anat"
    
    #    try to make output path and fail for reasons other than it already exists
    try:
        os.makedirs(out_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

        
    if not os.path.isfile(out_path + out_file + "+orig.HEAD"):

        #    need dicom path from _run_map_anat_
        dcm_path = subj.get_dcm_path(study,"anat","MPRAGE1")
        if not dcm_path:
            raise utils.CmdProcError("Dicom path not found")
            
        utils.prettyOut("Create anat volume...")

        #    dimon will dump stuff where we are running so change directory to output path
        os.chdir(out_path)

        if not os.path.isfile(out_path + "__proc.anat.tmp.to3d+orig.HEAD"):
            
            cmd = ( "Dimon -infile_pattern '" + dcm_path + "/*' -dicom_org -save_details dimon.details.txt -save_errors -gert_quit_on_err -gert_create_dataset -gert_to3d_prefix __proc.anat.tmp.to3d" )
            utils.Cmd(cmd).echo()
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)
    
            utils.purge_paths(out_path + "GERT_Reco_*")
                
            
        #    crop this to help with align_epi_anat
        if os.path.isfile(out_path + "__proc.anat.tmp.to3d+orig.HEAD"):
            cmd = ( "3dautobox -input " + out_path + "__proc.anat.tmp.to3d+orig -prefix " + out_path + "__proc.anat.tmp.box" )
            utils.Cmd(cmd).echo()
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)

                
        
        exit("REMOVE THIS STEP? OR leave for overall improvement of foreign data? Look again in saggital plane. Maybe put flag in function call")
        #    do some intensity correction to improve the skullstripping and align_epi_anat
        #    this output becomes our new base volume
        if os.path.isfile(out_path + "__proc.anat.tmp.box+orig.HEAD"):
            cmd = ("3dUniformize -anat " + out_path + "__proc.anat.tmp.box+orig -prefix " + out_path + out_file)
            utils.Cmd(cmd).echo()
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)        
        
            utils.purge_paths(out_path + "__proc.anat.tmp*", out_path + "p*.1D", out_path + "vt*.1D")
            
    
    
    

def skullstrip(subj,study,orig_vol=False):
    
    """
    skullstrip T1
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    out_path = study.path + "/" + subj.ursi.full + subj.visit_path + "/anat/"
    in_file  = os.path.join(out_path, study.label + "." + subj.ursi.short + ".anat+orig")
    o_file_unif = os.path.join(out_path, study.label + "." + subj.ursi.short + ".anat.unif+orig")
    o_file_sksp_tmp = os.path.join(out_path, study.label + "." + subj.ursi.short + ".anat.unif.sksp")
    o_file_sksp = os.path.join(out_path, study.label + "." + subj.ursi.short + ".anat.sksp")
    
    
    if not os.path.isfile(o_file_sksp + "+orig.HEAD"):

        utils.prettyOut("SkullStrip anat...")

        # unifize
        if not os.path.isfile(o_file_unif + ".HEAD"):
            utils.Cmd( "3dUnifize -GM -input " + in_file + " -prefix " + o_file_unif ).run()
        
        # skullstrip
        if not os.path.isfile(o_file_sksp_tmp + "+orig.HEAD"):
            utils.Cmd( "3dSkullStrip -input "+o_file_unif+" -prefix "+o_file_sksp_tmp+" -debug 1 -ld 23 -niter 777 -shrink_fac_bot_lim 0.75 -exp_frac 0.0666 -orig_vol -pushout -use_edge -touchup").run()
        
        # smooth and fill the sksp a bit
        if not os.path.isfile(o_file_sksp_tmp + ".mask+orig.HEAD"):
            utils.Cmd("3dmask_tool -input "+o_file_sksp_tmp+"+orig -dilate_input 2 -2 -prefix "+o_file_sksp_tmp + ".mask").run()

        # mask our original or unifized anat for a final sksp volume
        if orig_vol:
            utils.Cmd("3dcalc -a "+o_file_sksp_tmp + ".mask+orig.HEAD -b "+in_file+" -expr 'step(a)*b' -prefix " + o_file_sksp).run()
        else:
            utils.Cmd("3dcalc -a "+o_file_sksp_tmp + ".mask+orig.HEAD -b "+o_file_unif+" -expr 'step(a)*b' -prefix " + o_file_sksp).run()

        # clean up
        utils.purge_paths(o_file_unif+"*", o_file_sksp_tmp+"*", out_path+"*.anat.sksp_radrat.1D.dset")
            
            
        
        
def skullstripOLD(subj,study,orig_vol=False):
    
    """
    skullstrip T1
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    out_path = study.path + "/" + subj.ursi.full + subj.visit_path + "/anat/"
    in_file  = study.label + "." + subj.ursi.short + ".anat+orig"
    out_file = study.label + "." + subj.ursi.short + ".anat.sksp"
    
    
    if not os.path.isfile(out_path + out_file + "+orig.HEAD"):

        utils.prettyOut("SkullStrip anat...")
        #cmd = ( "3dSkullStrip -input " + out_path + in_file + " -prefix " + out_path + out_file + " -shrink_fac_bot_lim 0.8 -shrink_fac 0.8 -niter 300 ")
        
        #    unifize - don't do. It has very strange effect on T1W
        #cmd = ( "3dUnifize -overwrite -input  " + out_path + in_file + " -prefix " + out_path + "proc.anat.tmp.unifize")
        #utils.cmd_echo(cmd)
        #return_code = subprocess.check_call(cmd, shell=True)
        #if return_code:
            #exit(return_code)
        
        if orig_vol:
            cmd = ( "3dSkullStrip -orig_vol -input " + out_path + in_file + " -prefix " + out_path + out_file + " -shrink_fac_bot_lim 0.75 -shrink_fac 0.85 -niter 500 ")
            #cmd = ( "3dSkullStrip -orig_vol -input " + out_path + "proc.anat.tmp.unifize+orig -prefix " + out_path + out_file + " -shrink_fac_bot_lim 0.75 -shrink_fac 0.85 -niter 500 ")
        else:
            cmd = ( "3dSkullStrip -input " + out_path + in_file + " -prefix " + out_path + out_file + " -shrink_fac_bot_lim 0.75 -shrink_fac 0.85 -niter 500 ")
        utils.Cmd(cmd).echo()
        
        return_code = subprocess.check_call(cmd, shell=True)
        if return_code:
            exit(return_code)

        #utils.purge_paths(out_path + "proc.anat.tmp.unifize+orig*.1D")
            
            
            
    


def affine_align_to_Template(subj,study):
    
    """
    affine align T1 to Template (TT_N27)
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """        
    
    out_path = study.path + "/" + subj.ursi.full + subj.visit_path + "/anat/"
    in_file  = study.label + "." + subj.ursi.short + ".anat.sksp"
    
    
    if not os.path.isfile(out_path + in_file + "+tlrc.HEAD"):

        utils.prettyOut("affine_align_T1_to_Template anat...")
        
        #    sad limitation of this auto_tlrc
        os.chdir(out_path)
        cmd = ( "@auto_tlrc -no_ss -suffix 'NONE' -rmode quintic -base TT_N27+tlrc. -input " + in_file + "+orig -init_xform AUTO_CENTER" )
        utils.Cmd(cmd).echo()
        return_code = subprocess.check_call(cmd, shell=True)
        if return_code:
            exit(return_code)
            
    
    #    "-no_pre" flag was not cleaning up so it was removed and the "pre..." file are deleted manually here
    utils.purge_paths(out_path + "pre." + in_file + "*")
    
    
    
    


def NL_align_to_Template(subj,study):
    
    """
    Qwarp T1 to Template (TT_N27)
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    #Qwarp registration to the TT_N27+tlrc brain should be done by blurring both the base (TT_N27) and the source datasets.
    #-blur 3 3 -workhard:0:1

    import glob

    #    deposit junk in /tmp
    os.chdir("/tmp")
    
    
    #    build subject_path
    subject_path     = study.path + "/" + subj.ursi.full + "/" + subj.visit

    
    #    Assumes affine-aligned skullstripped input. Stops if that isn't available
    anat_sksp_file        = subject_path + "/anat/" + study.label +"."+ subj.ursi.short+".anat.sksp+tlrc"
    anat_sksp_unif_file    = subject_path + "/anat/" + study.label +"."+ subj.ursi.short+".anat.sksp.unif+tlrc"
    aff_3x4_tx            = subject_path + "/reg_matrices/" + study.label +"."+ subj.ursi.short+".T1_2tlrc.aff3x4.1D"
    final_WARP_file        = subject_path + "/reg_matrices/" + study.label +"."+ subj.ursi.short+".T1_2tlrc.WARP"
    anat_qwarp_file        = subject_path + "/anat/" + study.label +"."+ subj.ursi.short+".anat.sksp.qwarp"
    #    if you change template file you'll need to redo affine alignment to make sure it used the same template
    template_file        = utils.Env.afni_path + "/TT_N27+tlrc"
    #template_file        = utils.Env.afni_path + "/TT_avg152T1+tlrc"
    


    #    final out exists?
    if not os.path.isfile(anat_qwarp_file+"+tlrc.HEAD"):
        
        #    no final output, how about the final transform?
        if not os.path.isfile(final_WARP_file+"+tlrc.HEAD"):
        
        
            #    check inputs
            if not os.path.isfile(anat_sksp_file+".HEAD"):
                exit("\nMissing Input " + anat_sksp_file + "\n")
                
            
            #    Unifize skullstripped input
            if not os.path.isfile(anat_sksp_unif_file+".HEAD"):
                utils.prettyOut("Unifize skullstripped input")
                cmd = ( "3dUnifize -GM -prefix " + anat_sksp_unif_file + " -input " + anat_sksp_file )
                return_code = subprocess.check_call(cmd, shell=True)
                if return_code:
                    exit(return_code)
                        

            #    expects auto_tlrc already run before here
                    
            #    perhaps reg_matrices not there
            if not os.path.isdir(subject_path + "/reg_matrices/"):
                utils.prettyOut("Create directory: " + subject_path + "/reg_matrices/")
                cmd = ( "mkdir -p " + subject_path + "/reg_matrices" )
                return_code = subprocess.check_call(cmd, shell=True)
                if return_code:
                    exit(return_code)
                    
                    
            #    look for 3x4 T1 to tlrc matrix, if not there, create         
            if not os.path.isfile(aff_3x4_tx):
                utils.prettyOut("Create 3x4 T1 to tlrc transform")
                cmd = ( "cat_matvec " + anat_sksp_file + "::WARPDRIVE_MATVEC_FOR_000000 >  " + aff_3x4_tx )
                print("Create 3x4 T1 to tlrc transform")
                return_code = subprocess.check_call(cmd, shell=True)
                if return_code:
                    exit(return_code)        
            
            #    To nonlinear warp
            if not os.path.isfile(subject_path + "/reg_matrices/" + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp_WARP1+tlrc.HEAD"):
                utils.prettyOut("Create nonlinear warp")
                
                #    let's not bombard hosts seem to slow way down
                os.environ['OMP_NUM_THREADS'] = '6'
                
                #cmd = ( "3dQwarp -nodset -prefix " + anat_qwarp_file + " -duplo -useweight -blur 3 3 -base " + template_file + " -source " + anat_sksp_unif_file )
                # 20160308 - change warp flags given new information
                #cmd = ( "3dQwarp -nodset -prefix " + anat_qwarp_file + " -useweight -blur 3 3 -workhard:0:1 -base " + template_file + " -source " + anat_sksp_unif_file )
                # 20180319
                #cmd = ( "3dQwarp -nodset -prefix " + anat_qwarp_file + " -useweight -duplo -base " + template_file + " -source " + anat_sksp_unif_file )
                # 20190218
                cmd = ( "3dQwarp -nodset -prefix " + anat_qwarp_file + " -useweight -pblur 0.09 0.09 -Qfinal -workhard:5:8 -base " + template_file + " -source " + anat_sksp_unif_file )
                
                
                return_code = subprocess.check_call(cmd, shell=True)
                if return_code:
                    exit(return_code)
                    
                #    save the aff to tlrc qwarp so it can be catenated for other purposes in the future
                #    could just use the final WARP produced in next step but if we save this one we can catenate other things later without the extra interpolation of the warp
                os.rename(anat_qwarp_file+"_WARP+tlrc.HEAD",subject_path + "/reg_matrices/" + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp_WARP1+tlrc.HEAD")
                os.rename(anat_qwarp_file+"_WARP+tlrc.BRIK.gz",subject_path + "/reg_matrices/" + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp_WARP1+tlrc.BRIK.gz")
            
                #    clear this even though it should only be local to this python shell
                if 'OMP_NUM_THREADS' in os.environ:
                    del os.environ['OMP_NUM_THREADS']            
            
            
            #    create final native T1 to NL tlrc WARP
            #    20141030 - jl - As of 22 Oct 2014 this is not longer necessary as this concatenation can be moved into 3dNwarpApply call
            utils.prettyOut("create final native T1 to NL tlrc WARP")
            cmd = ( "3dNwarpCat -prefix " + final_WARP_file + " -warp1 " + subject_path + "/reg_matrices/" + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp_WARP1+tlrc -warp2 " + aff_3x4_tx  )
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)
            
            
            
        #    apply final nl transform to original anat.sksp > anat.sksp.qwarp+tlrc
        utils.prettyOut("apply final nl transform to original anat.sksp > anat.sksp.qwarp+tlrc")
        #cmd = ( "3dNwarpApply -interp 'wsinc5' -prefix " + anat_qwarp_file + " -nwarp " + final_WARP_file + "+tlrc -source " + subject_path + "/anat/" + study.label +"."+ subj.ursi.short+".anat.sksp+orig -master " + final_WARP_file  + "+tlrc" )
        cmd = ( "3dNwarpApply -interp 'wsinc5' -prefix " + anat_qwarp_file + " -nwarp " + final_WARP_file + "+tlrc -source " + subject_path + "/anat/" + study.label +"."+ subj.ursi.short+".anat.sksp+orig -master " + template_file )
        return_code = subprocess.check_call(cmd, shell=True)
        if return_code:
            exit(return_code)

        #    clean up
        trash = glob.glob(anat_sksp_unif_file+".*")
        for f in trash:
            os.remove(f)        


        
        utils.prettyOut("Completed qwarp T1 to tlrc")
    
    
    
    
    
def qa_NL_align_to_Template(subj,study):


    """
    QA display of Qwarp T1 to Template (TT_N27)
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    import glob
    
    
    os.chdir("/tmp")
    
    cmd = ( "3dedge3 -overwrite -prefix /tmp/" + study.label + "." +  subj.ursi.short + ".anat.sksp.qwarp.edge -input " + subj.get_path(study) + "anat/" + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp+tlrc" )
    print(cmd)
    return_code = subprocess.check_call(cmd, shell=True)
    if return_code:
        exit(return_code)
    
    
    
    # scale the volume for better contrast
    cmd = ( "@ScaleVolume -input /tmp/" + study.label + "." +  subj.ursi.short + ".anat.sksp.qwarp.edge+tlrc -prefix /tmp/" + study.label + "." +  subj.ursi.short + ".anat.sksp.qwarp.edge.sc -perc_clip 2 98 -val_clip 6 249")
    print(cmd)
    return_code = subprocess.check_call(cmd, shell=True)
    if return_code:
        exit(return_code)
        

    os.environ["AFNI_DETACH"] = "NO"

    cmd = ('afni -noplugins -com "SWITCH_UNDERLAY TT_N27" -com "SWITCH_OVERLAY ' + study.label + '.' + subj.ursi.short + '.anat.sksp.qwarp.edge.sc" -com "SEE_OVERLAY +" -com "SET_THRESHOLD A.0" -com "SET_PBAR_ALL A.+99 1.0" -com "SET_SPM_XYZ 6 12 14" ')
    print(cmd)
    return_code = subprocess.check_call(cmd, shell=True)
    if return_code:
        exit(return_code)
    
    os.environ["AFNI_DETACH"] = "YES"
    
    
    
    trash = glob.glob("/tmp/" + study.label + "." +  subj.ursi.short + ".anat.*")
    #print(trash)
    for f in trash:
        os.remove(f)
    

    
    
    
    
def generate_spm_tissue_classes(subj,study):
    
    
    """
    Invoke SPM tissue classification on T1
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    import glob
    
    
    
    input_path = subj.get_path(study) + "anat/"
    
    
    #    check for spm output
    c_name = sorted(glob.glob(input_path + "c*anat.nii"))
    #    check for final afni files
    a_name = sorted(glob.glob(input_path + "*.mask.60+orig.HEAD"))


    if len(a_name) != 3:
        
        if len(c_name) != 3:
        
            if not os.path.isfile(input_path + study.label + "." + subj.ursi.short + ".anat.nii"):
                
                if not os.path.isfile(input_path + study.label + "." + subj.ursi.short + ".anat.nii.gz"):
                    
                    cmd = ("3dresample -orient ASR -inset " + input_path + study.label + "." + subj.ursi.short + ".anat+orig -prefix " + input_path + study.label + "." + subj.ursi.short + ".anat.nii.gz")
                    utils.Cmd(cmd).echo()
                    return_code = subprocess.check_call(cmd, shell=True)
                    if return_code:
                        exit(return_code)
                
                cmd = ("gunzip " + input_path + study.label + "." + subj.ursi.short + ".anat.nii.gz")
                utils.Cmd(cmd).echo()
                return_code = subprocess.check_call(cmd, shell=True)
                if return_code:
                    exit(return_code)
                    

            
            if os.path.isfile(input_path + study.label + "." + subj.ursi.short + ".anat.nii"):
                
                utils.prettyOut("ANAT: generate spm tissue classes")
                
                #    read study settings. if these parameters are not set we won't send them to spm_new_segment()
                #    if they are, pass them all
                spm_settings = study.get_settings('spm_spatial_biasfwhm','spm_spatial_biasreg','spm_spatial_warp_cleanup')
                
                if (spm_settings['spm_spatial_biasfwhm'] != None):
                    cmd = ("matlab -nodisplay -nosplash -nojvm -singleCompThread -r \"spm_new_segment('" + input_path + study.label + "." + subj.ursi.short + ".anat.nii'," + spm_settings['spm_spatial_biasfwhm'] + "," + spm_settings['spm_spatial_biasreg'] + "," + spm_settings['spm_spatial_warp_cleanup'] + ")\"")
                else:
                    cmd = ("matlab -nodisplay -nosplash -nojvm -singleCompThread -r \"spm_new_segment('" + input_path + study.label + "." + subj.ursi.short + ".anat.nii')\"")
        
                utils.Cmd(cmd).echo()
                return_code = subprocess.check_call(cmd, shell=True)
                if return_code:
                    exit(return_code)
                    
                try:
                    os.remove(input_path + study.label + "." + subj.ursi.short + ".anat.nii")
                except:
                    raise
        
        
        
        
        
                
def convert_spm_tissue_classes(subj,study):
    
    """
    Takes output of SPM tissue classification on T1 and converts to AFNI format.
    Also, creates binary tissue masks.
    
    Separate this function for studies that roll their own.
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    import glob
    
    input_path = subj.get_path(study) + "anat/"
    if not os.path.exists(input_path):
        print(input_path)
        exit("Anat path not found")
    

    #    do three tissue types
    tissue_dict = dict(gm='1',wm='2',csf='3')
    for key in tissue_dict:
        
        out_file         = subj.get_path(study) + "anat/"+ study.label + "." + subj.ursi.short + ".anat.spm." + key
        anat_master_file    = subj.get_path(study) + "anat/"+ study.label + "." + subj.ursi.short + ".anat+orig"
        
        
    
        # already have what we need?
        if not os.path.isfile(out_file+".mask.60+orig.HEAD"):
            
            
            #    before we get started let's make sure we have anat_master_file
            if not os.path.exists(anat_master_file+".HEAD"):
                print(anat_master_file)
                exit("Anat file not found")
        
        
            utils.prettyOut("ANAT: convert spm tissue classes")
            print(("tissue:"+tissue_dict[key]))
            print(out_file)
            
            # don't change this order of operations
            if not os.path.isfile(out_file+".prob+tlrc.HEAD"):
                    
                f_name = glob.glob(input_path + "c"+tissue_dict[key]+"*")
                print("working file:")
                print(f_name)
                
                if len(f_name) != 1:
                    print(f_name)
                    exit("Ambiguous tissue file")

                cmd = ("3dresample -orient ASR -inset " + f_name[0] + " -prefix " + out_file + ".prob")
                return_code = subprocess.check_call(cmd, shell=True)
                if return_code:
                    exit(return_code)
                
            
                
            # change orient as needed and adjust origin to match sksp
            cmd = ("3drefit -view 'orig' -space 'ORIG' " + out_file + ".prob+tlrc")
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)
            
            
            cmd = ( "3drefit -atrcopy "+anat_master_file+" 'IJK_TO_DICOM' -atrcopy "+anat_master_file+" 'IJK_TO_DICOM_REAL' " + out_file +".prob+orig" )
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)
            
            cmd = ( "3drefit -duporigin " + anat_master_file + " " + out_file + ".prob+orig" )
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)
                
                            
                
            # create a binary of the mask at .60 threshold
            # another way to do this and include the partial volumeing would be to assign voxels to 1
            # if they are more than a third of probability and the others are less than a third
            # With the threshold method you will get overlap between tissue assignation                        
            cmd = ( "3dcalc -a " + out_file + ".prob+orig -prefix " + out_file + ".mask.60 -datum short -expr 'astep(a,0.60)'" )
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)
                
                
                
                
        del_list = glob.glob(input_path + "c"+tissue_dict[key]+"*")
        for f in del_list:
            os.remove(f)

            
        if os.path.isfile(subj.get_path(study) + "anat/m"+ study.label + "." + subj.ursi.short + ".anat.nii"):
            os.remove(subj.get_path(study) + "anat/m"+ study.label + "." + subj.ursi.short + ".anat.nii")
                
                
                
                
    
def pull_AA_spm_tissue_classes(subj,study):
    
    """
    Pulls SPM tissue classification from Auto_Analysis output.
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    import glob
    
    #    get session on run name from anat map
    dcm_path = subj.get_dcm_path(study,"anat","MPRAGE1")
    if not dcm_path:
        raise utils.CmdProcError("Dicom path not found")
            
    #run_name = dcm_path.split("/")[-1]
    session  = dcm_path.split("/")[-2]
    
    input_path = study.get_auto_analysis_path() + "triotim/" + subj.ursi.full + "/" + session + "/analysis/vbm/"
    if not os.path.exists(input_path):
        print(input_path)
        exit("Tissue class path not found")
    
    
    #    do three tissue types
    tissue_dict = dict(gm='1',wm='2',csf='3')
    for key in tissue_dict:
        
        out_file         = subj.get_path(study) + "anat/"+ study.label + "." + subj.ursi.short + ".anat.spm." + key
        anat_master_file    = subj.get_path(study) + "anat/"+ study.label + "." + subj.ursi.short + ".anat+orig"
        
        
        # already have what we need?
        if not os.path.isfile(out_file+".mask.60+orig.HEAD"):
            
            # don't change this order of operations
            if not os.path.isfile(out_file+".prob+tlrc.HEAD"):
                    
                f_name = glob.glob(input_path + "c"+tissue_dict[key]+"*")
                
                if len(f_name) != 1:
                    print(f_name)
                    exit("Ambiguous tissue file")

                cmd = ("3dresample -orient ASR -inset " + f_name[0] + " -prefix " + out_file + ".prob")
                print(cmd)
                return_code = subprocess.check_call(cmd, shell=True)
                if return_code:
                    exit(return_code)
                
            
                
            # change orient as needed and adjust origin to match sksp
            cmd = ("3drefit -view 'orig' -space 'ORIG' " + out_file + ".prob+tlrc")
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)
            
            
            cmd = ( "3drefit -atrcopy "+anat_master_file+" 'IJK_TO_DICOM' -atrcopy "+anat_master_file+" 'IJK_TO_DICOM_REAL' " + out_file +".prob+orig" )
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)        
            
            cmd = ( "3drefit -duporigin " + anat_master_file + " " + out_file + ".prob+orig" )
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)    
                
                            
                
            # create a binary of the mask at .60 threshold
            # another way to do this and include the partial volumeing would be to assign voxels to 1
            # if they are more than a third of probability and the others are less than a third
            # With the threshold method you will get overlap between tissue assignation                        
            cmd = ( "3dcalc -a " + out_file + ".prob+orig -prefix " + out_file + ".mask.60 -datum short -expr 'astep(a,0.60)'" )
            print(cmd)
            return_code = subprocess.check_call(cmd, shell=True)
            if return_code:
                exit(return_code)
                
                
            # some clean up
            del_list = glob.glob(out_file + "+orig*")
            for f in del_list:
                os.remove(f)

                
                
                
                
def qa_spm_tissue_classes(subj,study):

    """
    QA display of SPM tissue classification final masks.
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    # afni driver doesn't appear to respond to SET_SESSION, go to directory in OS
    os.chdir(os.path.join(subj.get_path(study),"anat") )
        
    # let afni know not to detach
    os.environ["AFNI_DETACH"] = "NO"
    
    # run the viewing command
    utils.Cmd('afni -noplugins \
    -com "OPEN_WINDOW A.coronalimage geom=320x320+889+457" \
    -com "SET_XHAIRS axialimage.OFF" \
    -com "SET_DICOM_XYZ 12 0 29" \
    -com "SWITCH_UNDERLAY '+ study.label + '.' + subj.ursi.short + '.anat" \
    -com "SWITCH_OVERLAY '+ study.label + '.' + subj.ursi.short + '.anat.spm.csf.mask.60" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 3000" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 300" \
    -com "SWITCH_OVERLAY '+ study.label + '.' + subj.ursi.short + '.anat.spm.wm.mask.60" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 3000" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 300" \
    -com "SWITCH_OVERLAY '+ study.label + '.' + subj.ursi.short + '.anat.spm.gm.mask.60" \
    -com "SEE_OVERLAY +"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY -"  \
    -com "SLEEP 500" \
    -com "SEE_OVERLAY +" ').run()

    
    os.environ["AFNI_DETACH"] = "YES"
    
    
    



    
def fs_recon_all(subj,study,flag_mprage=True):
    
    """
    Runs FreeSurfer recon-all on T1.
    This currently only uses raw dicom input.
    
    :param subj: utils.Subj
    :param study: utils.Study    
    :param flag_mprage: Should recon-all run with "-mprage" flag?
    
    """    
    
    try:
    
        # local definition of freesurfer path
        fs_path = os.path.join(study.path,"freesurfer",subj.visit)
        
        
        # freesurfer looks in environment for this
        os.environ['SUBJECTS_DIR'] = fs_path
        
        # check for fs output first
        if os.path.exists(os.path.join(fs_path,subj.ursi.full)):
            
            # check for partial completion
            if not os.path.isfile(os.path.join(fs_path,subj.ursi.full,"mri","brainmask.mgz")) and not os.path.isfile(os.path.join(fs_path,subj.ursi.full,"mri","T1.mgz")):
                exit("Half-baked Freesurfer output?")
            
            else:
                utils.prettyOut("Freesurfer recon_all output already exists")
                
        else:
            
            #    need dicom path from _run_map_anat_
            dcm_path = subj.get_dcm_path(study,"anat","MPRAGE1")
            if not dcm_path:
                raise utils.CmdProcError("Dicom path not found")
                
            dcm_one  = glob.glob( os.path.join(dcm_path,"SerieMR-*-0001*.dcm") )
            
            if len(dcm_one):
                
                # run recon_all
                cmd = utils.Cmd("recon-all -all -parallel -threads 9 -3T -i "+dcm_one[0]+" -subject "+subj.ursi.full)
                # 
                if flag_mprage:
                    cmd.append("-mprage")
                    
                cmd.run()
            
            else:
                exit("Failed to locate dicom file at " + dcm_path)

                        
    except:
        
        raise
    
    
    
    
def fs_clean_aseg(subj,study):
    
    """
    Apply a post-hoc clean to aseg mask. 
    This is appropriate for FreeSurfer version 5.2 but may be automatic in newer versions.
    
    :param subj: utils.Subj
    :param study: utils.Study    
    
    """    
    
    import shutil
    
    try:
    
        # local definition of freesurfer path
        fs_path = os.path.join(study.path,"freesurfer",subj.visit)
        
                
        # freesurfer looks in environment for this
        os.environ['SUBJECTS_DIR'] = fs_path

        # look for input first        
        if os.path.isfile(os.path.join(fs_path,subj.ursi.full,"mri","aseg.mgz")):
            
            #         
            utils.Cmd("apas2aseg --s "+subj.ursi.full).run()
            
            shutil.move( os.path.join(fs_path,subj.ursi.full,"mri","apas-aseg.mgz"), os.path.join(fs_path,subj.ursi.full,"mri","aseg.mgz") )
            

    except:
        
        raise
    
    
    
    
def fs_xhemi(subj,study):
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




    
def extract_raw_data(subj,study,run_label):

    """
    Extracts image for a run

    :param subj: utils.Subj
    :param study: utils.Study
    :param run_label: Run type with digit (e.g. 'FLAIR1'). Should be mapped in _run_map_anat_.

    """

    #run_settings = study.get_run_label_settings(run_label)
    out_path = subj.get_path(study) + run_label + "/"
    out_file_prefix = study.label + "." + subj.ursi.short + "." + run_label
    dcm_path = subj.get_dcm_path(study,"anat",run_label)

    # fail quietly if no dcm files found. Probably better to get run map instead and complain.
    if dcm_path:
    
        #    try to make output path and fail for reasons other than it already exists
        try:
            os.makedirs(out_path,0o770)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        if not os.path.isfile(out_path + out_file_prefix + ".nii.gz"):
            utils.Cmd( "dcm2niix -f "+ out_path + out_file_prefix + " -b y -g i -o " + out_path + " -v y -x y -z y" + dcm_path + "/*.dcm" ).run()
            
            # rename dcm2nii output
            #os.rename(out_path + subj.ursi.full + ".nii.gz", out_path + out_file_prefix + ".nii.gz")





class longitudinal:

    class mindboggle:
            
        def ants_cortical_thickness(subj, study, template_name='tpl.v1_v2'):
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
            import shutil
            
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
                        
