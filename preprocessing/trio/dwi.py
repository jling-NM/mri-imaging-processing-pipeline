"""@package dwi
Preprocessing module for diffusion imaging
Author: jling@mrn.org
DWI,DTI,DKI,etc.
"""

from mayerlab import utils   #
import os,errno,glob         #


# default output directory
proc_dir = 'dwi'


def run_qc_edc(subj,study):

    """
    QA display of EDDY results for inspection.

    :param subj: utils.Subj
    :param study: utils.Study

    """

    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : EDDY check" )
    
    os.chdir(subj.path_append(study, proc_dir))
    
    in_underlay =  subj.id_fn(study,"dwi.concat.nii.gz")
    if not os.path.isfile(in_underlay):
        in_underlay = subj.id_fn(study,"dwi.concat.nii")
        if not os.path.isfile(in_underlay):
            in_underlay = subj.id_fn(study,"dwi.concat.bc.nii.gz")
            if not os.path.isfile(in_underlay):
                raise utils.CmdProcError("Missing input:\n"+in_underlay)

    in_overlay = subj.id_fn(study,"dwi.concat.edc.nii.gz")
    if not os.path.isfile(in_overlay):
        in_overlay = subj.id_fn(study,"dwi.concat.bc.edc.nii.gz")
        if not os.path.isfile(in_overlay):
            raise utils.CmdProcError("Missing input:\n"+in_overlay)


    cmd = utils.Cmd('afni -com "SWITCH_UNDERLAY ' + in_underlay + '" -com "SWITCH_OVERLAY '+ in_overlay +'" -com "SEE_OVERLAY +" -com "SET_THRESHOLD A.0" -com "SET_PBAR_ALL A.+99 1.0" -com "SET_SPM_XYZ 0 0 0" -com "SEE_OVERLAY -" -com "SLEEP 500" -com "SEE_OVERLAY +" -com "SLEEP 500" -com "SEE_OVERLAY -" -com "SLEEP 500" -com "SEE_OVERLAY +"')
    cmd.run()
    
    
    
    
def get_edc_image_outlier_summary(subj,study):
    """
    Return outlier summary from eddy_outlier_map, and qc
    """
    dwi_path = subj.get_path(study) + "dwi"
    base_name   = study.label + "." + subj.ursi.short + "." + "dwi.concat.edc"
    
    import numpy as np
    
    # load outlier matrix
    try:
        
        i_file = os.path.join(dwi_path, base_name + '.eddy_outlier_map')
        if os.path.isfile(i_file):
            outlier_map = np.loadtxt( os.path.join(dwi_path, base_name + '.eddy_outlier_map'), skiprows=1, dtype='int')
            total_num_slices = outlier_map.shape[0] * outlier_map.shape[1]
    
            sum = np.sum(outlier_map, axis=1)
            pct = np.round((sum/outlier_map.shape[1])*100)
            r = np.array([list(range(outlier_map.shape[0])), sum, pct.astype(dtype='int')])
            
            outlier_slice_sum =np.sum(r[1])
            outlier_image_count1 = np.count_nonzero(r[1])
            outlier_image_count3 = np.count_nonzero(r[1]>3)
            outlier_slice_pct = (outlier_slice_sum/total_num_slices)*100
            outlier_image_count1_pct = (outlier_image_count1/outlier_map.shape[0])*100
            outlier_image_count3_pct = (outlier_image_count3/outlier_map.shape[0])*100
            
            return outlier_slice_sum, outlier_slice_pct, outlier_image_count1, outlier_image_count1_pct, outlier_image_count3, outlier_image_count3_pct, total_num_slices
            
        else:
            print("{} couldn't be located.".format(i_file))
            
    except:
        raise
        
        

def view_edc_image_outlier_summary(subj,study):
    """
    Generate outlier summary from eddy_outlier_map
    """
    # load outlier matrix
    try:
        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Gradient QC")
        
        ( outlier_slice_sum, outlier_slice_pct, outlier_image_count1, outlier_image_count1_pct, outlier_image_count3, outlier_image_count3_pct, total_num_slices ) = get_edc_image_outlier_summary(subj,study)
        print( "{} outlier artifacts found in {:0.2f}% of the slices.".format(outlier_slice_sum, outlier_slice_pct ) )
        print( "{} of the images({:0.2f}%) have 1 or more slice outliers. [Single Zipper].".format( outlier_image_count1, outlier_image_count1_pct ) )
        print( "{} of the images({:0.2f}%) have 3 or more slice outliers. [Multiple Zippers].".format( outlier_image_count3, outlier_image_count3_pct) )
            
    except:
        raise
    
    
    
    
    
    

def run_eddy_quad_qc(subj,study):
    """
    Generate eddy qc output
    """
    
    out_path    = subj.get_path(study) + "dwi"
    base_name   = study.label + "." + subj.ursi.short + "." + "dwi.concat.edc"
    
    try:
        if not os.path.isfile( os.path.join(out_path, 'qc','qc.json' )):
    
            utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Running QC")
    
            # go into dwi directory because of the way eddy_quad uses the full path in display
            os.chdir(out_path)
            utils.Cmd("eddy_quad " + base_name
                      + " -v --eddyIdx " + study.path + "/params/dwi_index.txt"
                      + " --eddyParams " + study.path + "/params/dwi_params.txt"
                      + " --mask " + base_name + ".mask.nii.gz"
                      + " --bvals " + study.label + "." + subj.ursi.short + "." + "dwi.concat.bvals"
                      + " --bvecs " + study.label + "." + subj.ursi.short + "." + "dwi.concat.edc.eddy_rotated_bvecs"
                      + " --output-dir qc " ).run()
    except:
        raise




def get_eddy_quad_qc_json(subj, study):
    """
    Returns stats from quad report in json format
    """
    
    try:
        out_path = subj.get_path(study) + "dwi"
        qc_file = os.path.join(out_path, 'qc','qc.json')
        if not os.path.isfile(qc_file):
            run_eddy_quad_qc(subj,study)
    
    
        # get stats
        import json

        with open(qc_file,'r') as f:
            qc_data = json.load(f)
            return qc_data

    except:
        raise


def view_eddy_quad_qc_stats(subj, study, raw=False):
    """
    Returns summary stats from quad report
    """
    
    try:
        qc_data = get_eddy_quad_qc_json(subj, study)
        
        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : EDDY QC Report Stats")
        print("Total Outliers (%): {:0.2f}".format(qc_data["qc_outliers_tot"]))
        print("Volume-to-Volume Motion:")
        print("   Average Abs motion (mm): {:0.2f}".format(qc_data["qc_mot_abs"]))
        print("   Average Rel motion (mm): {:0.2f}".format(qc_data["qc_mot_rel"]))
        
    except:
        raise







def generate_fsl_eddy_slspec(i_path_json, o_path):
    """
    Writes spec file using dcm2niix json file
    """
    import json
    import numpy as np
    
    # read json file
    with open(i_path_json) as f:    
        data = json.load(f)
        # get slice timing value
        slice_timing = np.array(data['SliceTiming'], dtype=float)
        # get mb factor
        mb = np.array(data['MultibandAccelerationFactor'], dtype=int)
        
    
    # sort to get index order
    sindx = np.argsort(slice_timing, axis=0, kind='mergesort')
    #sortedslicetimes = np.array(slice_timing)[sindx]
    
    # get mb 
    # essentially looking for number of slices collected at the same time. 2 slices have the same slicetime, mb=2, etc.
    #d = np.diff(sortedslicetimes, n=1)
    #mb = int(len(sortedslicetimes)/(len(d[d != 0])+1) )

    
    # generate slspec matrix
    slspec = np.reshape(sindx, (mb, int(len(sindx)/mb)), order='F')

    # save to o_file
    np.savetxt(o_path, slspec.T, fmt='%.d')
    
    
    
    



def create_bias_field(subj,study,run_label="DTI1",bias_bc_run_cnt=88,bias_32ch_run_cnt=88):
    """
    Create bias field.

    this is not particular to DWI, just here.

    :param subj: utils.Subj
    :param study: utils.Study
    :param run_label: Run type with digit (e.g. 'DTI2'). Should be mapped in _run_map_anat_. Default: 'DTI1'

    My implementation of Arvind's Auto-Analysis Bash/Matlab scripts

    """
    import glob

    out_path = os.path.join(subj.get_path(study),"bias")
    session_path = subj.get_session_path(study,"anat",run_label)


    #   let's confirm shit exists before we get all excited

    #   biasbc exists?
    #   note: bias run could have been repeated. Assume last is best and use that
    bias_bc_test = os.path.join(session_path,"bias_bc_*")
    bias_bc_paths = sorted(glob.glob(bias_bc_test))
    if len(bias_bc_paths) == 0:
        raise utils.CmdProcError("ERROR: Could not locate biasbc.\nCannot estimate bias field.")
    else:
        # we have a directory to use
        bias_bc_dcm_path = os.path.join(bias_bc_paths[-1],"*.dcm")
        # while the directory may be there, make sure the run is complete
        if len(glob.glob(bias_bc_dcm_path)) != bias_bc_run_cnt:
            raise utils.CmdProcError("ERROR: biasbc incomplete.\nCannot estimate bias field.")




    #    bias32 exists?
    bias_32ch_test = os.path.join(session_path,"bias_32ch_*")
    bias_32ch_paths = sorted(glob.glob(bias_32ch_test))
    if len(bias_32ch_paths) == 0:
        raise utils.CmdProcError("ERROR: Could not locate bias_32ch.\nCannot estimate bias field.")
    else:
        # we have a directory to use
        bias_32ch_dcm_path = os.path.join(bias_32ch_paths[-1],"*.dcm")
        # while the directory may be there, make sure the run is complete
        if len(glob.glob(bias_32ch_dcm_path)) != bias_32ch_run_cnt:
            raise utils.CmdProcError("ERROR: bias_32ch incomplete.\nCannot estimate bias field.")





    #    try to make output path and fail for reasons other than it already exists
    try:
        os.makedirs(out_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    if not os.path.isfile(os.path.join(out_path,"bias_full.nii.gz")):

        utils.prettyOut("DWI : Create bias correction image")

        # convert biasbc and bias32 to volume and dump in out_path
        if not os.path.isfile(os.path.join(out_path,"biasbc.nii.gz")):
            utils.Cmd( "dcm2niix -d N -e N -f N -g Y -n Y -o " + out_path + " -i N -p Y -r N -v Y -x N " + bias_bc_dcm_path ).run()


        if not os.path.isfile(os.path.join(out_path,"bias32ch.nii.gz")):
            utils.Cmd( "dcm2niix -d N -e N -f N -g Y -n Y -o " + out_path + " -i N -p Y -r N -v Y -x N " + bias_32ch_dcm_path ).run()


        # generate mysterious sform matrix that moves between biasbc and bias32
        o_file = os.path.join(out_path,"bias_bc_32ch.mtx")
        if not os.path.isfile(o_file):
            utils.Cmd(os.path.join(utils.Env.shared_script_path,"dicom_regc.sh") + " " + os.path.join(out_path,"biasbc.nii.gz") + " " + os.path.join(out_path,"bias32ch.nii.gz") + " " + o_file ).run()


        # move biasbc into bias32 grid using sform transform
        o_file = os.path.join(out_path,"biasbc_REG_INIT")
        if not os.path.isfile(o_file + ".nii.gz"):
            utils.Cmd( "flirt -in " + os.path.join(out_path,"biasbc.nii.gz") + " -ref " + os.path.join(out_path,"bias32ch.nii.gz") + " -out " + o_file + " -applyxfm -init " + os.path.join(out_path,"bias_bc_32ch.mtx") ).run()


        # skullstrip bias32
        if not os.path.isfile(os.path.join(out_path,"abias32ch_brain.nii.gz")):
            utils.Cmd( "bet " + os.path.join(out_path,"bias32ch.nii.gz") + " " + os.path.join(out_path,"abias32ch_brain") + " -m -f 0.4" ).run()


        # mask biasbc with our bias32 mask
        if not os.path.isfile(os.path.join(out_path,"abiasbc_brain.nii.gz")):
            utils.Cmd( "fslmaths " + os.path.join(out_path,"biasbc_REG_INIT.nii.gz") + " -mas " + os.path.join(out_path,"abias32ch_brain_mask") + " " + os.path.join(out_path,"abiasbc_brain") ).run()


        # invert our brain mask
        if not os.path.isfile(os.path.join(out_path,"abias_brain_invmask.nii.gz")):
            utils.Cmd( "fslmaths " + os.path.join(out_path,"abias32ch_brain_mask.nii.gz") + " -sub 1 -abs " + os.path.join(out_path,"abias_brain_invmask") ).run()


        # scale unmasked bias32
        if not os.path.isfile(os.path.join(out_path,"bias32ch_scaled.nii.gz")):
            utils.Cmd( "fslmaths " + os.path.join(out_path,"bias32ch.nii.gz") + " -add 0.1 " + os.path.join(out_path,"bias32ch_scaled") ).run()


        # scale unmasked biasbc, overwriting
        if not os.path.isfile(os.path.join(out_path,"biasbc_scaled.nii.gz")):
            utils.Cmd( "fslmaths " + os.path.join(out_path,"biasbc_REG_INIT.nii.gz") + " -add 0.1 " + os.path.join(out_path,"biasbc_scaled") ).run()

        # create a bias ratio data mask
        if not os.path.isfile(os.path.join(out_path,"abias_corr.nii.gz")):
            utils.Cmd( "fslmaths " + os.path.join(out_path,"biasbc_scaled.nii.gz") + " -div " + os.path.join(out_path,"bias32ch_scaled.nii.gz") + " -mul " + os.path.join(out_path,"abias32ch_brain_mask") + " " + os.path.join(out_path,"abias_corr") + " -odt float").run()


        # create some biasfit image
        if os.path.isfile(os.path.join(out_path,"abias_corr.nii.gz")):

            try:

                import nibabel as nb
                import numpy as np


                #    path to bias_corr input
                Fcorr_path  = os.path.join(out_path,"abias_corr.nii.gz")

                #    path to template
                Fbasis_path = os.path.join(study.path,"params","template_bias_basisn.nii")

                #    load Xcorr image
                Xcorr = nb.load(Fcorr_path)

                #    load template
                Xbasis = nb.load(Fbasis_path)

                #    shorten references to dimensions
                nx = Xcorr.shape[0]
                ny = Xcorr.shape[1]
                nz = Xcorr.shape[2]
                nt = Xbasis.shape[3]

                #    get bias_corr data
                Xcorr_data = Xcorr.get_data()

                #    vectorize bias_corr
                Xcorr_v = np.reshape(Xcorr_data,nx*ny*nz,1)

                #    filter out negative values and change datatype
                Xcorra = np.array(Xcorr_v[Xcorr_v > 0],dtype=np.double)

                #    reduce spatial dimensions of Xbasis and filter out same datapoints as we did for bias_corr
                Xbasis_data = Xbasis.get_data()
                Xbasis_v = np.reshape( Xbasis_data, (nx*ny*nz,nt), order='F')
                Xbasisa = np.array(Xbasis_v[np.where( Xcorr_v > 0 ),:],dtype=np.double)

                #    calculate bias fit
                u, d, v = np.linalg.svd(Xbasisa, full_matrices=False, compute_uv=1)

                #    d made a full matrix from just diagonal values
                dm = np.diagflat(d)
                #    shorten some references
                ut = np.transpose(u[0,:,:])
                v = v[0,:,:]
                vt= v.transpose()

                #    multiply all the matrices
                xa = np.dot(vt,np.dot(np.linalg.inv(dm),np.dot(ut,Xcorra)))

                #    multiply separately against bias template
                #    reshape for or image
                z = np.reshape(np.dot(Xbasis_v,xa),(nx,ny,nz), order='F')

                # write the data to new image using original bias_corr header
                newImg = nb.Nifti1Image(z,Xcorr.get_affine(), Xcorr.get_header())
                nb.save(newImg,os.path.join(out_path,"biasfit.nii.gz"))

            except:
                raise


        #  calculate bias_full image
        if not os.path.isfile(os.path.join(out_path,"bias_full.nii.gz")):
            utils.Cmd( "fslmaths " + os.path.join(out_path,"biasfit.nii.gz") + " -mul " + os.path.join(out_path,"abias_brain_invmask.nii.gz") + " -add " + os.path.join(out_path,"abias_corr.nii.gz") + " -kernel gauss 3 -fmedian -s 4 -thr 0 " + os.path.join(out_path,"bias_full")).run()


        # 11) remove everything save our bias_full
        utils.purge_paths(os.path.join(out_path,"abias*.nii.gz"),os.path.join(out_path,"biasfit.nii.gz"),os.path.join(out_path,"biasbc*"),os.path.join(out_path,"bias32ch*"),os.path.join(out_path,"biasbc_REG_INIT*"))







def extract_raw_data(subj,study, run_label):

    """
    Extracts image, bvals, and bvecs from dicom files.

    Using dcm2nii so works for scanners that dcm2nii works well with.
    NOTE: checks for dwi.concat.nii.gz before extracting individual run files again.

    :param subj: utils.Subj
    :param study: utils.Study
    :param run_label: Run type with digit (e.g. 'DTI2'). Should be mapped in _run_map_anat_.

    """
    
    #    check for any concat output before rerunning
    if not os.path.isfile(os.path.join(subj.path_append(study, proc_dir, study.label + "." + subj.ursi.short + "." + "dwi.concat.nii.gz"))):

        run_settings = study.get_run_label_settings(run_label)
        out_path = subj.get_path(study) + run_label + "/"
        out_file_prefix = study.label + "." + subj.ursi.short + "." + run_label

        #    create home for dwi output, do this here so it has time to create directory before it is needed below
        dwi_path = subj.append_path(study, proc_dir)
        try:
            os.makedirs(dwi_path,0o770)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        
        
        #    try to make output path and fail for reasons other than it already exists
        try:
            os.makedirs(out_path,0o770)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
                
        if not os.path.isfile(out_path + out_file_prefix + ".bvecs"):

            import glob, shutil

            utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Extract bvecs and bvals : " + run_label )

            #    need dicom path from _run_map_anat_
            dcm_path = subj.get_dcm_path(study,"anat",run_label)

            if not dcm_path:
                raise utils.CmdProcError("Dicom path not found")
                
            #    let's do a file count just so we can provide a more informative error.
            dcm_iter = glob.iglob(dcm_path + "/*.dcm")
            dcm_iter_len = len(list(dcm_iter))

            if(int(dcm_iter_len) != int(run_settings['time_images']) ):
                utils.errMsg( "Not enough dicom files at {}\n Expected {} but found {}".format(dcm_path + "/*.dcm",run_settings['time_images'],dcm_iter_len) )
                raise utils.CmdProcError("Not enough dicom files")

            if not os.path.isfile(out_path + out_file_prefix + ".nii.gz"):
                utils.Cmd( "dcm2niix -b y -g i -f "+out_file_prefix+" -z y -o " + out_path + " " + dcm_path ).run()

            if not os.path.isfile(out_path + out_file_prefix + ".bvec"):
                raise utils.CmdProcError("No bvec output from dcm2nii")

            #     prep bvec by transposing so they can be more easily manipulated
            if not os.path.isfile(out_path + out_file_prefix + ".bvecs"):
                utils.Cmd( "1d_tool.py -infile " + out_path + out_file_prefix + ".bvec -transpose -write " + out_path + out_file_prefix + ".bvecs" ).run()
                utils.purge_paths(out_path + out_file_prefix + ".bvec")

            #     rename bvals
            if not os.path.isfile(out_path + out_file_prefix + ".bvals"):
                os.rename(out_path + out_file_prefix + ".bval", out_path + out_file_prefix + ".bvals")

            bids_file = out_path + out_file_prefix + ".json"
            if os.path.isfile(bids_file):
                try:
                    shutil.copyfile( bids_file, os.path.join(dwi_path, "dwi"+run_label[-1]+".json") )
                except:
                    raise
                
            

def concatenate_runs(subj,study):

    """
    Concatenates multiple,expected DWI runs

    Checks that each run has been preprocessed before proceeding.

    :param subj: utils.Subj
    :param study: utils.Study

    """

    #    check for any concat output before rerunning
    if not os.path.isfile(os.path.join(subj.path_append(study, proc_dir, study.label + "." + subj.ursi.short + "." + "dwi.concat.nii.gz"))):

        import glob

        #    create home for dwi output
        out_path = subj.append_path(study, proc_dir)
        try:
            os.makedirs(out_path,0o770)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise



        #    check run settings and glob runs to make sure all have been preprocessed
        run_settings     = study.get_run_type_settings("DTI")
        comp_list        = sorted(glob.glob(subj.get_path(study) + "DTI?/" + study.label + "." + subj.ursi.short + "." + "DTI?.nii.gz"))


        #    if so, concatenate runs
        if ( len(comp_list) == int(run_settings['run_count']) ):

            utils.prettyOut("DWI : Concatentate Run Data")
            #    images; confirm sort order
            utils.Cmd( "3dTcat -prefix " + out_path + "/" + study.label + "." + subj.ursi.short + "." + "dwi.concat.nii " + subj.get_path(study) + "DTI?/" + study.label + "." + subj.ursi.short + "." + "DTI?.nii.gz" ).run()

            #    bvals
            #    concatenate and make sure no line feed in between
            utils.Cmd( "cat " + subj.get_path(study) + "DTI?/" + study.label + "." + subj.ursi.short + "." + "DTI?.bvals | awk 1 ORS=' ' > " + out_path + "/" + study.label + "." + subj.ursi.short + "." + "dwi.concat.bvals" ).run()

            #    bvecs
            utils.Cmd( "cat " + subj.get_path(study) + "DTI?/" + study.label + "." + subj.ursi.short + "." + "DTI?.bvecs  > " + out_path + "/" + study.label + "." + subj.ursi.short + "." + "dwi.concat.bvecs" ).run()

            # purge individual run directories to save space.
            if os.path.isfile(os.path.join(subj.path_append(study, proc_dir),study.label + "." + subj.ursi.short + "." + "dwi.concat.nii.gz")):
                utils.purge_paths(subj.get_path(study) + "DTI?")


        else:
            raise utils.CmdProcError("Cannot concatenate DTI runs before all individual runs have been preprocessed")




def verify_phase_encoding(subj,study):

    """

    :param subj: utils.Subj
    :param study: utils.Study

    """

    import json

    run_settings     = study.get_run_type_settings("DTI")

    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Phase Encoding Verification" )
    
    for run_num in range(1,int(run_settings['run_count'])+1):
        
        print("DWI"+str(run_num))
        
        i_path = os.path.join(subj.path_append(study, proc_dir),"dwi"+str(run_num)+".json")
        
        if os.path.isfile(i_path):
            with open(i_path) as sc_fh:
                sc = json.load(sc_fh)
                
                if( (sc['ProtocolName'].split('_')[3] == 'ap') & (sc['PhaseEncodingDirection'] != 'j-') ):
                    print("PhaseEncodingDirection Error: A>P encoding not implemented.")
                elif ( (sc['ProtocolName'].split('_')[3] == 'pa') & (sc['PhaseEncodingDirection'] != 'j') ):
                    print("PhaseEncodingDirection Error: P>A encoding not implemented.")
                else:
                    print("PhaseEncodingDirection Correct")
        else:
            print("No BIDS sidecar file.")
            
            






def round_bvals_file(subj,study):

    """
    Rounds b-values so they can be used with packages that expect that.
    Function uses subject and study input to find current bvals file.
    From that it creates a file in the same place and name except ".rounded"
    is added to end of filename.


    :param subj: utils.Subj
    :param study: utils.Study
    :returns path to rounded bvals file

    """
    
    in_bvals_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bvals"))
    if not os.path.isfile(in_bvals_path):
        raise utils.CmdProcError("Missing input:\n"+in_bvals_path)

    out_bvals_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bvals.rounded"))

    if not os.path.isfile(out_bvals_path):

        print(out_bvals_path)
        out_bvals_fh= open(out_bvals_path, "w")


        with open(in_bvals_path) as f:
            for line in f:
                for bval in line.strip().split(' '):
                    bval_i = int(round(float(bval+".0")/100)*100)
                    out_bvals_fh.write(str(bval_i)+" ")


        out_bvals_fh.close()

    return out_bvals_path



def get_bvecs_by_bvals(subj, study, bvalue):

    """
    Returns space-delimited bvec rows given a b-value

    :param subj: utils.Subj
    :param study: utils.Study

    """

    # we need rounded bvals file. Just call this function. It will check if it needs to create the file
    round_bvals_file(subj,study)


    in_bvals_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bvals.rounded") )
    if not os.path.isfile(in_bvals_path):
        raise utils.CmdProcError("Missing input:\n"+in_bvals_path)


    in_bvecs_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bvecs.adj"))
    if not os.path.isfile(in_bvecs_path):
        in_bvecs_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bvecs"))
        if not os.path.isfile(in_bvecs_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvecs_path)

    #cmd = utils.Cmd("cat " + in_bvals_path + " | tr ' ' '\n' | paste -d' ' - "+in_bvecs_path+" | grep '^"+bvalue+"' | cut -d' ' -f2-4")
    import subprocess

    cmd_txt = "cat " + in_bvals_path + " | tr ' ' '\\n' | paste -d' ' - " + in_bvecs_path + " | grep '^"+str(bvalue)+"' | cut -d' ' -f2-4"
    process = subprocess.Popen(cmd_txt,stdout=subprocess.PIPE, shell=True)
    cmd_result = process.stdout.read()

    return cmd_result





def eddy_cuda(subj, study, imain=None, mask=None, acqp=None, index=None, bvecs=None, bvals=None, topup=None, out=None, json=None, outlier_repol=True):
    """

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    study : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    from mayerlab import utils
    
    
    os.environ['SINGULARITYENV_CUDA_VISIBLE_DEVICES'] = utils.Gpu().get_available(25.0)
    if outlier_repol:
        utils.Cmd("/export/apps/linux-x86_64/singularity/bin/singularity exec --nv --cleanenv --bind /export:/export /export/research/analysis/human/amayer/shared/apps/docker/eddy_cuda.sif "
                  + "eddy_cuda "
                  + "--imain="+imain+" "
                  + "--mask="+mask+" "
                  + "--acqp="+acqp+" "
                  + "--index="+index+" "
                  + "--bvecs="+bvecs+" "
                  + "--bvals="+bvals+" "
                  + "--topup="+topup+" "
                  + "--out="+out+" "
                  + "--verbose "
                  + "--flm=quadratic "
                  + "--slm=none "
                  + "--fwhm=0 "
                  + "--niter=5 "
                  + "--interp=spline "
                  + "--resamp=jac "
                  + "--nvoxhp=1000 "
                  + "--ff=10 "
                  + "--cnr_maps "
                  + "--repol "
                  + "--ol_nstd=4 "
                  + "--ol_nvox=250 "
                  + "--ol_type=both "
                  + "--mporder=8 "
                  + "--niter=8 "
                  + "--fwhm=10,6,4,2,0,0,0,0 "
                  + "--s2v_niter=8 "
                  + "--s2v_lambda=1 "
                  + "--s2v_interp=trilinear "
                  + "--json="+json+" "
                  + "--estimate_move_by_susceptibility "
                  + "--mbs_niter=10 "
                  + "--mbs_lambda=10 "
                  + "--mbs_ksp=10").run()

    else:
        utils.Cmd("/export/apps/linux-x86_64/singularity/bin/singularity exec --nv --cleanenv --bind /export:/export /export/research/analysis/human/amayer/shared/apps/docker/eddy_cuda.sif "
                  + "eddy_cuda "
                  + "--imain="+imain+" "
                  + "--mask="+mask+" "
                  + "--acqp="+acqp+" "
                  + "--index="+index+" "
                  + "--bvecs="+bvecs+" "
                  + "--bvals="+bvals+" "
                  + "--topup="+topup+" "
                  + "--out="+out+" "
                  + "--verbose "
                  + "--flm=quadratic "
                  + "--slm=none "
                  + "--fwhm=0 "
                  + "--niter=5 "
                  + "--interp=spline "
                  + "--resamp=jac "
                  + "--nvoxhp=1000 "
                  + "--ff=10 "
                  + "--cnr_maps "
                  + "--mporder=8 "
                  + "--niter=8 "
                  + "--fwhm=10,6,4,2,0,0,0,0 "
                  + "--s2v_niter=8 "
                  + "--s2v_lambda=1 "
                  + "--s2v_interp=trilinear "
                  + "--json="+json+" "
                  + "--estimate_move_by_susceptibility "
                  + "--mbs_niter=10 "
                  + "--mbs_lambda=10 "
                  + "--mbs_ksp=10").run()
        


def eddy_openmp(subj, study, imain=None, mask=None, acqp=None, index=None, bvecs=None, bvals=None, topup=None, out=None, json=None, outlier_repol=True):
    """

    Parameters
    ----------
    subj : TYPE
        DESCRIPTION.
    study : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # how useful to include? Will only work with basic eddy functionality
            #    clear this even though it should only be local to this python shell
            # if 'OMP_NUM_THREADS' in os.environ:
            #     del os.environ['OMP_NUM_THREADS']



    
            # #    specifying defaults explicitly to avoid future surprises
            # #    i confirmed that the output of this matches implicit defaults at this time
            # #     call this "float" and it will later be converted
            # # eddy_openmp
            # utils.Cmd( "eddy_openmp "
            #     + "--imain=" + out_path + out_file_prefix + ".nii.gz "
            #     + "--mask="  + out_path + out_file_prefix + ".edc.mask.nii.gz "
            #     + "--acqp="  + study.path + "/params/dwi_params.txt "
            #     + "--index=" + study.path + "/params/dwi_index.txt "
            #     + "--bvecs=" + out_path + study.label + "." + subj.ursi.short + "." + "dwi.concat.bvecs "
            #     + "--bvals=" + out_path + study.label + "." + subj.ursi.short + "." + "dwi.concat.bvals "
            #     + "--topup=" + out_path + "topup_results "
            #     + "--out="   + out_path + out_file_prefix + ".edc "
            #     + "--verbose "
            #     + "--niter=5 "
            #     + "--fwhm=0 "
            #     + "--flm=quadratic "
            #     + "--resamp=jac ").run()





def susceptibility_eddy_correction(subj, study, b0_indices = ['0', '44']):

    """
    This corrects susceptibility-induced distortions with FSL EDDY
    using reversed phase-encoding polarities

    use the first and second run, just hardcoded since it is a one-time application by default but allow other input to accomadate bad acquisitions
    
    :param subj: utils.Subj
    :param study: utils.Study

    """
    import glob

    out_path          = subj.path_append(study, proc_dir)
    out_file_prefix   = study.label + "." + subj.ursi.short + "." + "dwi.concat"

    if not (len(glob.glob(os.path.join(out_path, out_file_prefix + ".edc.nii*" )))):

        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Eddy current and motion correction" )

        # verify params input first
        file_path_verify = [study.path + '/params/dwi_params.txt',study.path + '/params/dwi_index.txt']
        for fp in file_path_verify:
            if not os.path.isfile(fp):
                raise utils.CmdProcError("Missing Param input: " + fp)


        #    create b0 volume for topup if not already there
        b0_filepath = os.path.join(out_path, "topup.b0.nii.gz")
        extract_topup_b0(subj, study, b0_filepath, b0_indices)

        
        #    run topup for susceptibility estimate
        #    WARNING: Eddy cannot handle periods in the results file names
        #    TIME: ~ 5 mins on 4 core
        if not os.path.isfile(os.path.join(out_path, "topup_results_fieldcoef.nii.gz") ):

            utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Running Topup")

            utils.Cmd( "topup --verbose"
                    + " --imain=" + b0_filepath
                    + " --datain=" + study.path + "/params/dwi_params_b0.txt"
                    + " --config=b02b0_4.cnf"
                    + " --out=" + os.path.join(out_path, "topup_results")
                    + " --iout=" + os.path.join(out_path, "topup.dewarped_b0")
                    + " --logout=" + os.path.join(out_path, "topup_log") ).run()

        

        #    create dwi mask
        if not os.path.isfile( os.path.join(out_path, out_file_prefix + ".edc.mask.nii.gz" )):

            # may want to change this to use either only the first b0 or all the b0 with Tmean; RMMV
            # compare results
            utils.prettyOut("Create dwi mask")

            #    average susceptibility distortion corrected b0 volumes
            if not os.path.isfile(  os.path.join(out_path, "topup.dewarped_b0.Tmean.nii.gz" )):
                utils.prettyOut("Create b0 Tmean")
                utils.Cmd( "fslmaths " +  os.path.join(out_path, "topup.dewarped_b0") + " -Tmean " + os.path.join(out_path, "topup.dewarped_b0.Tmean") + " -odt short").run()

            
            #    mask the average susceptibility distortion corrected b0 mean volume
            # fsl mask, not very good
            #cmd = ( "bet " + out_path + out_file_prefix + ".topup.dewarped_b0.Tmean.nii.gz "  + out_path + out_file_prefix + ".topup.dewarped_b0 -f 0.30 -g 0 -n -m" )
            #    afni mask, better
            utils.Cmd("3dautomask -dilate 1 -prefix " + os.path.join(out_path, out_file_prefix + ".edc.mask.nii.gz") + " " + os.path.join(out_path, "topup.dewarped_b0.Tmean.nii.gz") ).run()


        

        #   run Eddy/Topup for distortion and motion
        #   DWI only support jac? eddy: msg=These data do not support least-squares resampling
        if not os.path.isfile( os.path.join(out_path, out_file_prefix + ".edc.nii.gz" )):

            utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Eddy/Topup for distortion and motion...Smoke, if ya got 'em")

            imain=os.path.join(out_path, out_file_prefix + ".nii.gz")
            mask=os.path.join(out_path, out_file_prefix + ".edc.mask.nii.gz")
            acqp=os.path.join(study.path,"params/dwi_params_b0.txt")
            index=os.path.join(study.path, "params/dwi_index_b0.txt")
            bvecs=os.path.join(out_path, study.label + "." + subj.ursi.short + "." + "dwi.concat.bvecs")
            bvals=os.path.join(out_path, study.label + "." + subj.ursi.short + "." + "dwi.concat.bvals")
            topup=os.path.join(out_path, "topup_results")
            out=os.path.join(out_path, out_file_prefix + ".edc")
            json=os.path.join(out_path, "dwi1.json")
                
            # dwi2 will be without outlier replacement
            eddy_cuda(subj, study, imain=imain, mask=mask, acqp=acqp, index=index, bvecs=bvecs, bvals=bvals, topup=topup, out=out, json=json, outlier_repol=True)


            # convert float datatype that comes out of EDDY
            if os.path.isfile( os.path.join(out_path, out_file_prefix + ".edc.nii.gz" )):
                utils.Cmd("fslmaths " + os.path.join(out_path, out_file_prefix + ".edc.nii.gz") + " " + os.path.join(out_path, out_file_prefix + ".edc") + " -odt short").run()



        #    transpose and rename EDDY rotated bvecs for afni programs
        #i_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.eddy_rotated_bvecs")
        i_bvecs_path = subj.path_append(study, proc_dir) + out_file_prefix + ".edc.eddy_rotated_bvecs"
        o_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvecs.adj")
        if not os.path.isfile( o_bvecs_path ):

            try:

                utils.prettyOut("Transpose\n" + i_bvecs_path + "\nto\n" + o_bvecs_path)

                from numpy import loadtxt,savetxt
                # transpose and rename at the same time
                savetxt(o_bvecs_path,loadtxt(i_bvecs_path).transpose(),fmt='%.6f')

            except:
                raise


        #    clean up
        utils.purge_paths(out_path + out_file_prefix + ".edc.float.nii.gz")

        exit('dwi.py')








def get_eddy_movement_data(subj,study,use_bias_in=False, out_units='radians'):

    """
    extracts and cleans dwi motion parameters output from FSL EDDY
    
    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#HelpOnRotBvec
    Understanding eddy output:
        eddy_output.eddy_parameters : This is a text file with one row for each volume and one column for each parameter. 
                                      The first six columns correspond to subject movement starting with three translations followed by three rotations. 
                                      The remaining columns pertain to the EC-induced fields and the number and interpretation of them will depend of which EC model was specified.


        eddy_output.eddy_movement_rms : A summary of the "total movement" in each volume is created by calculating the displacement of each voxel and then averaging the squares of those displacements across all intracerebral voxels (as determined by --mask and finally taking the square root of that. The file has two columns where the first contains the RMS movement relative the first volume and the second column the RMS relative the previous volume.
    
    
    
    Put output in a format similar to 3dvolreg to avoid confusion.
    i.e.
    roll pitch yaw dS dL dP

                    roll  = rotation about the I-S axis }
                    pitch = rotation about the R-L axis } degrees CCW
                    yaw   = rotation about the A-P axis }
                      dS  = displacement in the Superior direction  }
                      dL  = displacement in the Left direction      } mm
                      dP  = displacement in the Posterior direction }

    ASSUMING RPI
        rot_z rot_x rot_y trans_z trans_x trans_y
        5     3     4     2       0       1


    R = Y X Z where Z is rotation about the I-S axis, X is rotation about the R-L axis, and Y is rotation about the A-P axis

    this does not give the same output layout as 3dvolreg but puts "edc" in file name
    Rotations in this file are given in radians unlike 3dvolreg which is in degrees

    :param subj: utils.Subj
    :param study: utils.Study

    """

    #  define where movement file live for this study and visit
    out_path = study.movement_path + "/" + subj.visit

    #  make sure path exists or make it
    try:
        os.makedirs(out_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Extract motion parameters" )
    #    if not file, complain
    dwi_path = subj.get_path(study) + "dwi/"

    if use_bias_in:
        in_file_path = dwi_path + study.label + "." + subj.ursi.short + "." + "dwi.concat.bc.edc.eddy_parameters"
    else:
        in_file_path = dwi_path + study.label + "." + subj.ursi.short + "." + "dwi.concat.edc.eddy_parameters"



    if not os.path.isfile(in_file_path):
        raise utils.CmdProcError("Missing " + in_file_path)

    else:
        import numpy as np

        #    read rigid body alignment columns from eddy data file
        dwi_mov_data = np.loadtxt(in_file_path)[:,0:6]
        #    i don't like how reference image line is sometimes really small numbers so i force them to zero
        dwi_mov_data[0,:] = [0,0,0,0,0,0]
        
        if out_units == 'degrees':
            # convert radians to degrees
            dwi_mov_data[:,3:] = utils.Convert.Motion.radiansToDegrees(dwi_mov_data[:,3:])
            
        #    reorder; see notes above
        #    write out; place "edc" in file name to denote that this is not 3dvolreg output. e.g. rotations are in radians
        np.savetxt(out_path + "/" + subj.ursi.short + ".DWI.edc.Movement.Regressor.1D", dwi_mov_data[:,(5,3,4,2,0,1)], delimiter=' ', fmt='%1.6f')






def imReg_3d(subj,study):

    """
    This is the old school eddy current correction
    You might use this is you have multiple/high b-value data
    but no distortion images for EDDY and TOPUP

    If you have single, low b-value data you are better off
    with 3dvolreg which can be found in preprocessing_12ch.py
    scripts

    :param subj: utils.Subj
    :param study: utils.Study

    """

    import glob

    out_path        = subj.get_path(study) + "dwi/"
    out_file_prefix = study.label + "." + subj.ursi.short + "." + "dwi.concat"


    #    fsl eddy_correct and align
    if not (len(glob.glob(out_path + out_file_prefix + ".ecc.nii*" ))):

        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Eddy current and motion correction")


        if not os.path.isfile( out_path + out_file_prefix + ".ecc.nii.gz" ):
            utils.Cmd( "eddy_correct " + out_path + out_file_prefix + ".nii.gz " + out_path + out_file_prefix + ".ecc 0" ).run()


        utils.Cmd( "1d_tool.py -overwrite -infile " + out_path + out_file_prefix + ".bvecs -transpose -write " + out_path + out_file_prefix + ".bvecs.trans" ).run()

        utils.Cmd( "fdt_rotate_bvecs " + out_path + out_file_prefix + ".bvecs.trans " + out_path + out_file_prefix + ".bvecs.trans.adj " + out_path + out_file_prefix + ".ecc.ecclog" ).run()

        utils.Cmd( "1d_tool.py -overwrite -infile " + out_path + out_file_prefix + ".bvecs.trans.adj -transpose -write " + out_path + out_file_prefix + ".bvecs.adj" ).run()


    if not os.path.isfile(out_path + out_file_prefix + ".ecc.mask.nii.gz"):
        utils.Cmd("3dautomask -prefix " + out_path + out_file_prefix + ".ecc.mask.nii " + out_path + out_file_prefix + ".ecc.nii.gz").run()


    ##    create a masked version of the dewarped concatenated dwi
    #    also scaling data to short to save space, confirm that post-processing code does not complain about this
    if not os.path.isfile(out_path + out_file_prefix + ".ecc.masked.nii.gz"):

        utils.prettyOut("Mask distortion corrected data")

        cmd = utils.Cmd("fslmaths " + out_path + out_file_prefix + ".ecc "
                       + "-mas " + out_path + out_file_prefix + ".ecc.mask "
                       + out_path + out_file_prefix + ".ecc.masked "
                       + "-odt short ")

        cmd.run()


#        if not (len(glob.glob(out_path + out_file_prefix + ".ecc.masked.nii*" ))):
#            #    mask the end results
#            cmd = ( "bet " + out_path + out_file_prefix + ".ecc " + out_path + out_file_prefix + ".ecc.masked -f 0.4 -g 0 -n -m" )
#            utils.cmd_run(cmd)
#
        #fslmaths MMA.99373.dwi.concat.ecc.nii -mas MMA.99373.dwi.concat.ecc_brain_mask.nii.gz MMA.99373.dwi.concat.ecc.masked


    #    clean up
    utils.purge_paths(out_path + out_file_prefix + ".bvecs.trans*", out_path + out_file_prefix + "+orig*" )







def generate_dki(subj,study):

    """
    Diffusion Kurtosis Imaging

    :param subj: utils.Subj
    :param study: utils.Study

    """


    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DKI" )


    out_path        = subj.get_path(study) + "dwi/dki"
    try:
        os.makedirs(out_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


    # i'll input masked data along with the mask because the unmasked data is a large file
    in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.masked.nii.gz")
    if not os.path.isfile(in_dwi_path):
        in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.masked.nii.gz")
        if not os.path.isfile(in_dwi_path):
            in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.ecc.masked.nii.gz")
            if not os.path.isfile(in_dwi_path):
                raise utils.CmdProcError("Missing input:\n"+in_dwi_path)


    in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.mask.nii.gz")
    if not os.path.isfile(in_dwi_mask_path):
        in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.mask.nii.gz")
        if not os.path.isfile(in_dwi_mask_path):
            in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.ecc.mask.nii.gz")
            if not os.path.isfile(in_dwi_mask_path):
                raise utils.CmdProcError("Missing input:\n"+in_dwi_mask_path)


    in_bvals_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvals")
    if not os.path.isfile(in_bvals_path):
        raise utils.CmdProcError("Missing input:\n"+in_bvals_path)


    in_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvecs.adj")
    if not os.path.isfile(in_bvecs_path):
        in_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvecs")
        if not os.path.isfile(in_bvecs_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvecs_path)



    #
    o_file = os.path.join(out_path,"tmp.dwi.nii")

    if not os.path.isfile( o_file ):
        import shutil
        shutil.copyfile(in_dwi_path,o_file+'.gz')
        utils.Cmd( "gzip -d " + o_file+'.gz').run()

    #    run dki script
    #     Rscript /export/research/analysis/human/amayer/shared/apps/R/fit_dki.R --file_path=/export/research/analysis/human/amayer/tbi_07272/MMA_STUDY/dwi.analysis/dki/data --bvecs=bvecs.dat --bvals=bvals.dat --dwi=data.nii
    utils.Cmd( "Rscript --no-save --no-restore /export/research/analysis/human/amayer/shared/apps/R/fit_dki.R --file_path=" + out_path + " --bvecs=" + in_bvecs_path + " --bvals=" + in_bvals_path + " --dwi=" + o_file).run()

    #     clean up
    utils.purge_paths(o_file)






def denoise_dwi_nlmeans(subj,study):

    """
    Denoises DWI images using Non-Local Means (NLMEANS)

    :param subj: utils.Subj
    :param study: utils.Study

    """

    from os.path import join
    import nibabel as nib
    from dipy.denoise.nlmeans import nlmeans
    from dipy.denoise import noise_estimate


    data_path = subj.get_path(study) + "dwi"
    in_file   = join(data_path, study.label + "." + subj.ursi.short + "." + "dwi.concat")

    #    load image
    img = nib.load(in_file)
    aff = img.get_affine()
    #    pull data from image object
    data = img.get_data()

    #    noise estimate per image ( Standard deviation estimation from local patches http://nipy.org/dipy/reference/dipy.denoise.html )
    est_sigma = noise_estimate.estimate_sigma(data, disable_background_masking=False)
    #    calculate denoised dwi
    den = nlmeans(data, sigma=est_sigma, mask=None, patch_radius=1, block_radius=5, rician=True)
    #    write out denoised volume
    nib.save(nib.Nifti1Image(den, aff), 'denoised.rician.est_sigma.nii.gz')







#============================================================================
#
#
#============================================================================
def apply_bias_correction(subj,study):

    """
    Applies a bias correction image to DWI data (dwi.concat.bc.nii.gz).

    :param subj: utils.Subj
    :param study: utils.Study

    My implementation of Arvind's Auto-Analysis Bash/Matlab scripts

    """

    data_path = subj.get_path(study) + "dwi"
    data_file = study.label + "." + subj.ursi.short + "." + "dwi.concat.nii.gz"


    if not os.path.isfile( os.path.join(data_path,study.label + "." + subj.ursi.short + "." + "dwi.concat.bc.nii.gz") ):

        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Applying bias correction")

        #    extract b0
        out_file_path = os.path.join(data_path,"vol0000.nii.gz")
        if not os.path.isfile( out_file_path ):
            print("extract b0")
            utils.Cmd("fslroi " + os.path.join(data_path,data_file) + " " + out_file_path + " 0 1").run()


        #    vox grid align bias_full to dti
        out_file_path = os.path.join(data_path,"bias_dti_dic.mtx")
        if not os.path.isfile( out_file_path ):
            print("align bias_full to dti")
            utils.Cmd(os.path.join(utils.Env.shared_script_path,"dicom_regc.sh") + " " + os.path.join(subj.get_path(study),"bias","bias_full.nii.gz") + " " + os.path.join(data_path,"vol0000.nii.gz") +" "+ out_file_path ).run()


        out_file_path = os.path.join(data_path,"bias_dti_reg.nii.gz")
        if not os.path.isfile( out_file_path ):
            utils.Cmd("flirt -in " + os.path.join(subj.get_path(study),"bias","bias_full.nii.gz") + " -ref " + os.path.join(data_path,"vol0000.nii.gz") +" -out "+ out_file_path + " -applyxfm -init " + os.path.join(data_path,"bias_dti_dic.mtx") ).run()
            
        #    apply bias correction to dti
        out_file_path = os.path.join(data_path,study.label + "." + subj.ursi.short + "." + "dwi.concat.bc.nii.gz")
        if not os.path.isfile( out_file_path ):
            utils.Cmd("fslmaths " + os.path.join(data_path,data_file) + " -mul " + os.path.join(data_path,"bias_dti_reg.nii.gz") + " " + out_file_path + " -odt short ").run()


        #    clean up
        utils.purge_paths( os.path.join(data_path,"vol0000.nii.gz"), os.path.join(data_path,"bias_dti*") )







def spatial_normalization(subj,study,cost_func="lpc",big_move=True,cmass="-cmass cmass",opts="",purge=False,b0index="0"):

    """
    Performs affine registration and outputs transform dwi2tlrc.aff12.1D

    Expects masked b0 image will already be available.

    :param subj: utils.Subj
    :param study: utils.Study
    :param cost_func: valid cost functions for 3dAllineate
    :param big_move: [True|False]
    :param cmass:
    :param opts: Any other valid 3dAllineate flags to pass along
    :param purge: By default this function will not overwrite transforms as those may have been checked. Use this if you don't want to manually delete the previous transform.
    :param b0index: Brik index of b0 image to use for registration. Default = 0

    """


    # if alignment matrix exists and not told to blow it away
    if purge or not os.path.isfile(subj.get_transform_path(study) + subj.id_fn(study,"dwi2tlrc.aff12.1D")):


        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Spatial Normalization")

        # just so user understands
        if purge:
            utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : Purging previous transform matrix")
            utils.purge_paths(subj.get_transform_path(study) + "/*dwi*.1D")

        import shutil, glob

        input_fn_stem = "dwi.concat.bc.edc.masked"
        if not os.path.isfile(subj.path_append(study, proc_dir) + subj.id_fn(study,input_fn_stem+".nii.gz")):
            input_fn_stem = "dwi.concat.edc.masked"
            if not os.path.isfile(subj.path_append(study, proc_dir) + subj.id_fn(study,input_fn_stem+".nii.gz")):
                input_fn_stem = "dwi.concat.ecc.masked"
                if not os.path.isfile(subj.path_append(study, proc_dir) + subj.id_fn(study,input_fn_stem+".nii.gz")):
                    raise utils.CmdProcError("Missing input:\n"+input_fn_stem)



        # check for inputs
        input_file = os.path.join(subj.path_append(study, proc_dir),subj.id_fn(study, input_fn_stem+".nii.gz" ))

        if not os.path.isfile( input_file ):
            raise utils.CmdProcError("Cannot spatially normalize without input.\n"+input_file)


        #    create alignment workspace
        workspace_path = os.path.join(subj.path_append(study, proc_dir),"regTemp")
        try:
            os.makedirs(workspace_path,0o770)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        #    go in there cuz afni script is too dumb to use full paths correctly
        current_path = os.getcwd()
        os.chdir( workspace_path )


        ##    copying over working files
        for filename in glob.glob(os.path.join(subj.get_path(study,"anat"), subj.id_fn(study,"anat.sksp+*"))):
            shutil.copy(filename, workspace_path)

        #    here we copy only the b0 image of interest from the masked input
        cmd = utils.Cmd("fslroi "+os.path.join(subj.path_append(study, proc_dir),subj.id_fn(study,input_fn_stem+".nii.gz")) + " " +os.path.join(workspace_path,"b0") + " " + b0index + " 1")
        cmd.run()


        #    do alignment, sadly afni script not made callable from python
        if big_move:
            big_move_str = "-big_move"
        else:
            big_move_str = ""

        utils.Cmd("align_epi_anat.py \
                -epi2anat \
                -anat " + study.label + "." + subj.ursi.short + ".anat.sksp+orig \
                -epi " + os.path.join(workspace_path,"b0.nii.gz") + " \
                -epi_base 0 \
                -epi_strip None \
                -volreg off \
                -deoblique on \
                -anat_has_skull no \
                -tshift off \
                -cmass cmass "+big_move_str+" \
                -cost "+cost_func+" \
                -tlrc_apar " + subj.id_fn(study,"anat.sksp+tlrc") + " \
                -Allineate_opts '-source_automask+4'").run()

        # TODO: trap errors and clean up "regTemp" directory

        # create home for matrix if not there
        if not os.path.isdir(subj.get_transform_path(study)):
            try:
                os.makedirs(subj.get_transform_path(study),0o770)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        #    rename matrices
        utils.Cmd( "cat_matvec -ONELINE b0_al+orig.HEAD::ALLINEATE_MATVEC_B2S_000000 > " + subj.get_transform_path(study) + study.label + "." + subj.ursi.short + ".dwi2T1.aff12.1D" ).run()

        shutil.copy("b0_al_tlrc_mat.aff12.1D",subj.get_path(study,"tx") + subj.id_fn(study,"dwi2tlrc.aff12.1D"))
        # TODO: trap errors and clean up "regTemp" directory


        #    clean up
        os.chdir( current_path )
        shutil.rmtree(workspace_path, ignore_errors=True)






def qa_gradients(subj,study):

    """
    QA display of DWI for inspection.

    :param subj: utils.Subj
    :param study: utils.Study

    """

    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Check gradient images" )

    # check for inputs
    input_file_name =  subj.id_fn(study,"dwi.concat.nii.gz")
    if not os.path.isfile( os.path.join(subj.path_append(study, proc_dir),input_file_name) ):
        input_file_name =  subj.id_fn(study,"dwi.concat.nii")
        if not os.path.isfile( os.path.join(subj.path_append(study, proc_dir),input_file_name) ):
            input_file_name = subj.id_fn(study,"dwi.concat.bc.nii.gz")
            if not os.path.isfile( os.path.join(subj.path_append(study, proc_dir),input_file_name) ):
                raise utils.CmdProcError("Missing input:\n"+input_file_name)

    cmd = utils.Cmd('afni -com "OPEN_WINDOW axialimage geom=416x500+800+0 mont=3x4:3" -com "OPEN_WINDOW saggitalimage geom=500x500+100+650" \
         -com "OPEN_WINDOW coronalimage geom=500x300+800+650" -com "OPEN_WINDOW sagittalgraph" \
         -com "SWITCH_UNDERLAY ' + input_file_name + '" \
         -com "SET_XHAIRS OFF" ')

    # first go into that directory
    os.chdir(subj.path_append(study, proc_dir))

    cmd.run()







def qa_spatial_normalization(subj,study,b0index="0"):

    """
    QA display of DWI spatial normalization for inspection.

    This display DWI over subject's non-linear structural warp.

    :param subj: utils.Subj
    :param study: utils.Study
    :param b0index: Brik index of image to display.
    """

    import shutil, glob

    workspace_path = "/tmp"

    #    try doing this in machine tmp since this is all temporary output
    os.chdir(workspace_path)


    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : spatial normalization check" )

    if not os.path.isfile(subj.get_path(study,"tx") + subj.id_fn(study,"dwi2tlrc.aff12.1D")):
        raise utils.CmdProcError("Missing transform matrix. Perhaps you need to first spatially normalize?")


    if not os.path.isfile("VERIFY.dwi.tmp+tlrc.HEAD"):


        #
        in_file_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.masked.nii.gz")
        if not os.path.isfile(in_file_path):
            in_file_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.masked.nii.gz")
            if not os.path.isfile(in_file_path):
                in_file_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.ecc.masked.nii.gz")
                if not os.path.isfile(in_file_path):
                    raise utils.CmdProcError("Missing input:\n"+in_file_path)



        #    do SBRef due to nice contrast for visual check
        nl_warp_src        = subj.get_path(study,"tx") + study.label + "." + subj.ursi.short + ".anat.sksp.qwarp_WARP1+tlrc"
        t1_to_tlrc_mat    = subj.get_path(study,"tx") + study.label + "." + subj.ursi.short + ".T1_2tlrc.aff3x4.1D"
        epi_2_T1_aff_mat    = subj.get_path(study,"tx") + study.label + "." + subj.ursi.short + ".dwi2T1.aff12.1D"

        #    the reference volume output will be aligned to
        #    anat.sksp.qwarp+tlrc
        #    which is non-linear warp to TT_N27
        #    so reference volume output should look good over that
        #
        #    source is assuming that first brick is b0
        cmd = utils.Cmd("3dNwarpApply " +
                " -source " + in_file_path + "'["+b0index+"]'" +
                " -nwarp '" + nl_warp_src + " " + t1_to_tlrc_mat + " " + epi_2_T1_aff_mat + "'" +
                " -master " + subj.get_path(study,"anat") + study.label + "." + subj.ursi.short+".anat.sksp.qwarp+tlrc" +
                " -interp wsinc5" +
                " -prefix VERIFY.dwi.tmp" )

        cmd.run()

#        linear
#        cmd = ( "3dAllineate " +
#                " -master " + subj.get_path(study,"anat") + subj.id_fn(study,"anat.sksp+tlrc ") +
#                " -1Dmatrix_apply " + subj.get_path(study,"tx") + study.label + "." + subj.ursi.short + ".dwi2tlrc.aff12.1D" +
#                " -input " + in_file_path + "'[0]'" +
#                " -floatize" +
#                " -interp trilinear" +
#                " -final wsinc5"+
#                " -prefix VERIFY.dwi.tmp" )
#
#        utils.cmd_run(cmd)


    ##    copying over working files
    try:
        for filename in glob.glob(os.path.join(subj.get_path(study,"anat"), study.label + "." + subj.ursi.short + ".anat.sksp.qwarp+tlrc*")):
            shutil.copy(filename, workspace_path)
#
#
#
#        ##    copying over working files from anat
#        cmd = utils.Cmd("cp " + subj.get_path(study) + "anat/" + study.label+"."+subj.ursi.short+".anat.sksp.qwarp+tlrc* "  + workspace_path )
#        cmd.run()
#
    except:
        raise


    utils.Cmd( "3dedge3 -overwrite -prefix VERIFY.dwi.tmp.edge -input VERIFY.dwi.tmp+tlrc.HEAD" ).run()

    # scale the volume for better contrast
    utils.Cmd( "@ScaleVolume -input VERIFY.dwi.tmp.edge+tlrc -prefix VERIFY.dwi.tmp.edge.sc -perc_clip 2 98 -val_clip 6 249").run()

    #    Don't let afni session detach or it will close before user can see images
    os.environ["AFNI_DETACH"] = "NO"

    utils.Cmd('afni -com "SWITCH_UNDERLAY ' + study.label + '.' + subj.ursi.short + '.anat.sksp.qwarp" -com "SWITCH_OVERLAY VERIFY.dwi.tmp.edge.sc" -com "SEE_OVERLAY +" -com "SET_THRESHOLD A.0" -com "SET_PBAR_ALL A.+99 1.0" -com "SET_SPM_XYZ 6 12 14" ').run()

    #    return this setting to lab default
    os.environ["AFNI_DETACH"] = "YES"

    #    clean up
    utils.purge_paths("/tmp/" + study.label + "." +  subj.ursi.short + ".anat.*","/tmp/VERIFY.dwi.*")




def qa_dwi2T1(subj,study,b0index="0"):

    """
    QA display of DWI spatial normalization for inspection.

    This display DWI over subject's non-linear structural warp.

    :param subj: utils.Subj
    :param study: utils.Study
    :param b0index: Brik index of image to display.
    """

    import shutil, glob

    workspace_path = "/tmp"

    #    try doing this in machine tmp since this is all temporary output
    os.chdir(workspace_path)


    utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : to T1 check" )

    if not os.path.isfile(subj.get_path(study,"tx") + subj.id_fn(study,"dwi2T1.aff12.1D")):
        raise utils.CmdProcError("Missing transform matrix. Perhaps you need to first spatially normalize?")


    if not os.path.isfile("VERIFY.dwi.tmp+tlrc.HEAD"):

        #
        in_file_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.masked.nii.gz")
        if not os.path.isfile(in_file_path):
            in_file_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.masked.nii.gz")
            if not os.path.isfile(in_file_path):
                in_file_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.ecc.masked.nii.gz")
                if not os.path.isfile(in_file_path):
                    raise utils.CmdProcError("Missing input:\n"+in_file_path)



        utils.Cmd( "3dAllineate " +
                " -master " + subj.get_path(study,"anat") + subj.id_fn(study,"anat+orig ") +
                " -1Dmatrix_apply " + subj.get_path(study,"tx") + study.label + "." + subj.ursi.short + ".dwi2T1.aff12.1D" +
                " -input " + in_file_path + "'[0]'" +
                " -floatize" +
                " -interp trilinear" +
                " -final wsinc5"+
                " -prefix VERIFY.dwi.tmp" ).run()


    ##    copying over working files
    try:
        for filename in glob.glob(os.path.join(subj.get_path(study,"anat"), study.label + "." + subj.ursi.short + ".anat+orig*")):
            shutil.copy(filename, workspace_path)

    except:
        raise


    utils.Cmd( "3dedge3 -overwrite -prefix VERIFY.dwi.tmp.edge -input VERIFY.dwi.tmp+orig.HEAD" ).run()

    # scale the volume for better contrast
    utils.Cmd( "@ScaleVolume -input VERIFY.dwi.tmp.edge+orig -prefix VERIFY.dwi.tmp.edge.sc -perc_clip 2 98 -val_clip 6 249").run()

    #    Don't let afni session detach or it will close before user can see images
    os.environ["AFNI_DETACH"] = "NO"

    utils.Cmd('afni -com "SWITCH_UNDERLAY ' + study.label + '.' + subj.ursi.short + '.anat" -com "SWITCH_OVERLAY VERIFY.dwi.tmp.edge.sc" -com "SEE_OVERLAY +" -com "SET_THRESHOLD A.0" -com "SET_PBAR_ALL A.+99 1.0" -com "SET_SPM_XYZ 6 12 14" ').run()

    #    return this setting to lab default
    os.environ["AFNI_DETACH"] = "YES"

    #    clean up
    utils.purge_paths("/tmp/" + study.label + "." +  subj.ursi.short + ".anat.*","/tmp/VERIFY.dwi.*")














def generate_zhang_noddi(subj,study,use_bias_in=False):

    """
    Zhang implementation of NODDI

    :param subj: utils.Subj
    :param study: utils.Study

    """

    out_noddi_path = os.path.join(subj.path_append(study, proc_dir),"zhang_noddi")

    #    try to make noddi_path and fail for reasons other than it already exists
    try:
        os.makedirs(out_noddi_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    os.chdir(out_noddi_path)

    # zhang script doesn't run if there is a zero-length FittedParams.mat file from a crash or cancel. Clean up.
    fit_file_path = os.path.join(out_noddi_path,"FittedParams.mat")
    if (os.path.isfile(fit_file_path)) and (os.path.getsize(fit_file_path) == 0):
        os.remove(fit_file_path)
            
    
    if not os.path.isfile(os.path.join(out_noddi_path,"FittedParams.mat")):

        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : generate_zhang_noddi" )
            
        # i'll input masked data along with the mask because the unmasked data is a large file
        in_dwi_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.edc.nii.gz"))
        if not os.path.isfile(in_dwi_path):
            in_dwi_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bc.edc.nii.gz"))
            if not os.path.isfile(in_dwi_path):
                in_dwi_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.ecc.nii.gz"))
                if not os.path.isfile(in_dwi_path):
                    raise utils.CmdProcError("Missing input:\n"+in_dwi_path)

        
        in_dwi_mask_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.edc.mask.nii.gz"))
        if not os.path.isfile(in_dwi_mask_path):
            in_dwi_mask_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bc.edc.mask.nii.gz"))
            if not os.path.isfile(in_dwi_mask_path):
                in_dwi_mask_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.ecc.mask.nii.gz"))
                if not os.path.isfile(in_dwi_mask_path):
                    raise utils.CmdProcError("Missing input:\n"+in_dwi_mask_path)


        in_bvals_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bvals"))
        if not os.path.isfile(in_bvals_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvals_path)


        in_bvecs_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bvecs.adj"))
        if not os.path.isfile(in_bvecs_path):
            in_bvecs_path = subj.path_append(study, proc_dir, subj.id_fn(study,"dwi.concat.bvecs"))
            if not os.path.isfile(in_bvecs_path):
                raise utils.CmdProcError("Missing input:\n"+in_bvecs_path)



        # decompress dwi file into out_path
        import gzip
        # what files must be decompressed
        file_paths_gzip = [in_dwi_path, in_dwi_mask_path]

        for fp in file_paths_gzip:

            # first get filename without gz for output
            fp_filename = '.'.join(fp.split(os.sep)[-1].split('.')[0:-1])
            fp_path     = os.path.join(out_noddi_path,fp_filename)

            if not os.path.isfile(fp_path):

                utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : decompress input" )

                decompressedFile = gzip.GzipFile(fp, mode='rb')
                with open(fp_path, 'wb') as fp_path:
                    fp_path.write(decompressedFile.read())


        # what files must be copied
        file_paths_copy = [in_bvals_path, in_bvecs_path]
        import shutil
        for fp in file_paths_copy:
            fp_filename = fp.split(os.sep)[-1]
            dst_file_path = os.path.join(out_noddi_path,fp_filename)
            if not os.path.isfile(dst_file_path):
                shutil.copy(fp, dst_file_path)


        # unwrap above selected filenames for matlab script
        dwi_f      = ".".join(in_dwi_path.split(os.path.sep)[-1].split(".")[0:-1])
        dwi_mask_f = ".".join(in_dwi_mask_path.split(os.path.sep)[-1].split(".")[0:-1])
        bvals_f    = in_bvals_path.split(os.path.sep)[-1]
        bvecs_f    = in_bvecs_path.split(os.path.sep)[-1]


        #    call matlab function ( this should clean up )
        # dwi_noddi(noddi_path,ursi,study_label,dwi_input,dwi_mask,bvals_path,bvecs_path)
        cmd = utils.Cmd("matlab -nodisplay -nosplash -nojvm -r \"dwi_noddi('" + out_noddi_path + "','" + subj.ursi.short + "','" + study.label + "','" + dwi_f + "','" + dwi_mask_f + "','" + bvals_f + "','" + bvecs_f + "')\"")
        cmd.run()



    #  confirm matlab output before proceeding
    if not os.path.isfile(os.path.join(out_noddi_path,'FittedParams.mat')):
        raise utils.CmdProcError("Missing dwi_noddi output. NODDI incomplete.\n")


    # once everything is done running, re-indent this section and remove import gzip
    import gzip

    # do some output renaming to match other methodologies
    ren_dict = {'noddi_ficvf': 'fic', 'noddi_fiso': 'fiso', 'noddi_fmin': 'min', 'noddi_kappa': 'kappa', 'noddi_odi': 'odi'}
    for ren_key in ren_dict:

        n_path = os.path.join(out_noddi_path,subj.id_fn(study,ren_key+".nii"))
        z_path = os.path.join(out_noddi_path,subj.id_fn(study,ren_dict[ren_key]+".nii.gz"))

        if not os.path.isfile(z_path):

            if not os.path.isfile(n_path):
                raise utils.CmdProcError("Stopping on missing dwi_noddi output: "+n_path+"\nNODDI incomplete.\n")

            # zip it
            with open(n_path, 'rb') as nii_file:
                with gzip.open(z_path, 'wb') as zipped_file:
                    zipped_file.writelines(nii_file)

            # remove decompressed version
            if os.path.isfile(z_path):
                os.remove(n_path)



    # calculate some variants of the volume fractions

    # fic prime = (1-fiso)*fic
#    out_f = os.path.join(out_noddi_path,subj.id_fn(study,"fic_iso_adj.nii.gz"))
#    if not os.path.isfile(out_f):
#        cmd = utils.Cmd("3dcalc -a "+ os.path.join(out_noddi_path,subj.id_fn(study,"fiso.nii.gz")) +" -b "+ os.path.join(out_noddi_path,subj.id_fn(study,"fic.nii.gz")) +" -prefix "+out_f+" -expr '(1-a)*b'")
#        cmd.run()
#

    # basic fec calculation
    # fec       = step(odi)*(1-fic)
    # the step(odi) reduces map to voxel not zero in odi. Because we need some way to mask and fiso=0 would not be a good mask.
#    out_f = os.path.join(out_noddi_path,subj.id_fn(study,"fec.nii.gz"))
#    if not os.path.isfile(out_f):
#        cmd = utils.Cmd("3dcalc  -a "+ os.path.join(out_noddi_path,subj.id_fn(study,"fic.nii.gz")) +" -b "+ os.path.join(out_noddi_path,subj.id_fn(study,"odi.nii.gz")) +" -prefix "+out_f+" -expr '(step(b)*(1-a))'")
#        cmd.run()


    # fec       = (1-fiso)*(step(odi)*(1-fic))
    # the step(b) reduces map to voxel not zero in fic. Because we need some way to mask and fiso=0 would not be a good mask.
    # note: we are using unadjusted fic here
    #out_f = os.path.join(out_noddi_path,subj.id_fn(study,"fec.nii.gz"))
    #if not os.path.isfile(out_f):
    #    cmd = utils.Cmd("3dcalc -a "+ os.path.join(out_noddi_path,subj.id_fn(study,"fiso.nii.gz")) +" -b "+ os.path.join(out_noddi_path,subj.id_fn(study,"fic.nii.gz")) +" -c "+ os.path.join(out_noddi_path,subj.id_fn(study,"odi.nii.gz")) +" -prefix "+out_f+" -expr '(1-a)*(step(c)*(1-b))'")
    #    cmd.run()









def generate_amico_noddi(subj,study,protocol_path):

    """
    Daducci implementation of NODDI

    :param subj: utils.Subj
    :param study: utils.Study
    :param protocol_path: full path to data common and reused for noddi in this study ( e.g. kernals )

    """

    out_noddi_path = os.path.join(subj.path_append(study, proc_dir),"amico_noddi")

    #    try to make noddi_path and fail for reasons other than it already exists
    try:
        os.makedirs(out_noddi_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    os.chdir(out_noddi_path)




    if not os.path.isfile(os.path.join(out_noddi_path,subj.id_fn(study,"fiberDirs.nii.gz"))):


        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : generate_amico_noddi" )


        # i'll input masked data along with the mask because the unmasked data is a large file
        in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.masked.nii.gz")
        if not os.path.isfile(in_dwi_path):
            in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.masked.nii.gz")
            if not os.path.isfile(in_dwi_path):
                in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.ecc.masked.nii.gz")
                if not os.path.isfile(in_dwi_path):
                    raise utils.CmdProcError("Missing input:\n"+in_dwi_path)


        in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.mask.nii.gz")
        if not os.path.isfile(in_dwi_mask_path):
            in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.mask.nii.gz")
            if not os.path.isfile(in_dwi_mask_path):
                in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.ecc.mask.nii.gz")
                if not os.path.isfile(in_dwi_mask_path):
                    raise utils.CmdProcError("Missing input:\n"+in_dwi_mask_path)







        # decompress dwi file into out_path
        import gzip
        # what files must be decompressed
        file_paths_gzip = [in_dwi_path, in_dwi_mask_path]

        for fp in file_paths_gzip:

            # first get filename without gz for output
            fp_filename = '.'.join(fp.split(os.sep)[-1].split('.')[0:-1])
            fp_path     = os.path.join(out_noddi_path,fp_filename)

            if not os.path.isfile(fp_path):

                utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : decompress input" )

                decompressedFile = gzip.GzipFile(fp, mode='rb')
                with open(fp_path, 'wb') as fp_path:
                    fp_path.write(decompressedFile.read())



        # we need rounded bvals file for this metric. Just call this function. It will check if it needs to create the file
        round_bvals_file(subj,study)


        # scheme file
        generate_camino_scheme(subj,study)

        in_scheme_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.camino.scheme")
        if not os.path.isfile(in_scheme_path):
            raise utils.CmdProcError("Missing input:\n"+in_scheme_path)



        # copy scheme file into out_path
        import shutil,subprocess
        shutil.copy(in_scheme_path, os.path.join(out_noddi_path,subj.id_fn(study,"dwi.concat.camino.scheme")))



          # calculate amico noddi
        if not os.path.isfile(os.path.join(out_noddi_path,"FIT_dir.nii")):

            #    call matlab function
            # cmd = utils.Cmd("matlab -nodisplay -nosplash -nojvm -r \"dwi_amico('" + out_noddi_path + "','"+'.'.join(in_dwi_path.split(os.sep)[-1].split('.')[0:-1])+"','"+'.'.join(in_dwi_mask_path.split(os.sep)[-1].split('.')[0:-1])+"','MMA." + subj.ursi.short + ".dwi.concat.camino.scheme')\"")
            # ( protocol_path, data_in_path, dwiFilename, maskFilename, schemeFilename, data_out_path );
            cmd = utils.Cmd("matlab -nodisplay -nosplash -nojvm -r \"dwi_amico('" + protocol_path + "','" + out_noddi_path + "','"+'.'.join(in_dwi_path.split(os.sep)[-1].split('.')[0:-1])+"','"+'.'.join(in_dwi_mask_path.split(os.sep)[-1].split('.')[0:-1])+"','MMA." + subj.ursi.short + ".dwi.concat.camino.scheme','" + out_noddi_path + "')\"")
            cmd.run()

            # clean up duplicate inputs    first or use a more specific pattern
            utils.purge_paths(os.path.join(out_noddi_path,"MMA.*"))



        # do some output renaming
        dict_file_names_rename = { 'FIT_dir.nii':subj.id_fn(study,"fiberDirs.nii"),'FIT_ICVF.nii':subj.id_fn(study,"fic.nii"),'FIT_ISOVF.nii':subj.id_fn(study,"fiso.nii"),'FIT_OD.nii':subj.id_fn(study,"odi.nii") }
        for key in dict_file_names_rename:
            if os.path.isfile(key):
                # grab the provenence
                p = subprocess.Popen(['3dnotes',key], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (stdoutdata, stderrdata) = p.communicate()
                prov = stdoutdata.replace("\n----- HISTORY -----\n","")

                # copy the matlab output with 3dcopy otherwise AFNI will complain about the disconnect in the header
                # also deals with header ID conflict caused by matlab nifti write used in AMICO code
                cmd = utils.Cmd("3dcopy " + key + " " +dict_file_names_rename[key] )
                cmd.run()
                # update the header so we know how data was processed
                cmd = utils.Cmd('3dnotes -HH "'+prov+'" ' + dict_file_names_rename[key] + '.gz' )
                cmd.run()

                # out with the old
                os.remove(key)





        # first, fix upper limit on fiso volume to 1.0; values greater than that are probably in CSF
        cmd = utils.Cmd("3dcalc -overwrite -a " + subj.id_fn(study,"fiso.nii.gz") + " -prefix " + subj.id_fn(study,"fiso.nii.gz") + " -expr 'a+(astep(a,1)*(1-a))'" )
        cmd.run()




    # calculate some variants of the volume fractions

#    # fic prime = (1-fiso)*fic
#    out_f = subj.id_fn(study,"fic_iso_adj.nii.gz")
#    if not os.path.isfile(out_f):
#        cmd = utils.Cmd("3dcalc -a "+subj.id_fn(study,"fiso.nii.gz")+" -b "+subj.id_fn(study,"fic.nii.gz")+" -prefix "+out_f+" -expr '(1-a)*b'")
#        cmd.run()
#
#
#    # fec       = (1-fiso)*(step(odi)*(1-fic))
#    # the step(b) reduces map to voxel not zero in fic. Because we need some way to mask and fiso=0 would not be a good mask.
#    # note: we are using unadjusted fic here
#    out_f = subj.id_fn(study,"fec.nii.gz")
#    if not os.path.isfile(out_f):
#        cmd = utils.Cmd("3dcalc -a "+subj.id_fn(study,"fiso.nii.gz")+" -b "+subj.id_fn(study,"fic.nii.gz")+" -c "+subj.id_fn(study,"odi.nii.gz")+" -prefix "+out_f+" -expr '(1-a)*(step(c)*(1-b))'")
#        cmd.run()







def generate_nonlinear_DT(subj,study,use_bvals='all'):

    """
    Standard AFNI nonlinear diffusion tensor fit.

    :param subj: utils.Subj
    :param study: utils.Study
    :param use_bvals: which gradients to use [some b-value,'all']

    """

    in_bvals_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvals")
    if not os.path.isfile(in_bvals_path):
        raise utils.CmdProcError("Missing input:\n"+in_bvals_path)


    in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.masked.nii.gz")
    if not os.path.isfile(in_dwi_path):
        in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.masked.nii.gz")
        if not os.path.isfile(in_dwi_path):
            in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.ecc.masked.nii.gz")
            if not os.path.isfile(in_dwi_path):
                raise utils.CmdProcError("Missing input:\n"+in_dwi_path)


    in_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvecs.adj")
    if not os.path.isfile(in_bvecs_path):
        in_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvecs")
        if not os.path.isfile(in_bvecs_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvecs_path)



    # by default use all gradients but if some b-value list is passed, use just those
    if use_bvals == 'all':

        brik_selector_str = ''
        import shutil
        shutil.copy(in_bvecs_path,subj.path_append(study, proc_dir)+"3dDWItoDT.bvecs.use_bvals")

    else:

        select_bvals = [0]
        select_bvals.append(int(use_bvals))
        brik_selector_ind = []

        brik_indx = 0
        with open(in_bvals_path) as f:
            for line in f:
                for bval in line.strip().split(' '):

                    bval_i = int(round(float(bval+".0")/100)*100)
                    if bval_i in select_bvals:
                        brik_selector_ind.append(brik_indx)

                    brik_indx += 1


        brik_selector_str = "'["+ ','.join(map(str, brik_selector_ind)) + "]'"



        new_bvecs_file =  open(subj.path_append(study, proc_dir)+"3dDWItoDT.bvecs.use_bvals", "w")

        with open(in_bvecs_path) as f:
            brik_indx = 0
            for line in f:
                if brik_indx in brik_selector_ind:
                    new_bvecs_file.write(line)
                brik_indx += 1


        new_bvecs_file.close()

    out_dti_path = os.path.join(subj.path_append(study, proc_dir),"dti")
    #    try to make dti path and fail for reasons other than it already exists
    try:
        os.makedirs(out_dti_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


    out_dti_pre  = subj.id_fn(study,"DTI_NL.nii")
    #if not os.path.isfile(os.path.join(out_dti_path,subj.id_fn(study,"DTI_NL_FA.nii.gz"))):
    if not len( glob.glob(os.path.join(out_dti_path,subj.id_fn(study,"DTI_NL_FA.nii*"))) ):

        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Nonlinear Calculation of Tensor" )
        utils.prettyOut("Using "+use_bvals+" gradients for model")


        # create 3dDWItoDT version of bvecs by dropping initial b0 gradient
        utils.Cmd("tail -n +2 "+ subj.path_append(study, proc_dir)+"3dDWItoDT.bvecs.use_bvals"+" > " + subj.path_append(study, proc_dir)+"3dDWItoDT.bvecs" ).run()
        utils.Cmd("3dDWItoDT -prefix  " + os.path.join(out_dti_path,out_dti_pre) +
                            " -nonlinear -max_iter 10 -max_iter_rw 5 -reweight -eigs -sep_dsets -mean_b0" +
                            " " + subj.path_append(study, proc_dir)+"3dDWItoDT.bvecs " +
                            in_dwi_path + brik_selector_str).run()

    utils.purge_paths(subj.path_append(study, proc_dir)+"3dDWItoDT.bvecs*")






def generate_camino_scheme(subj,study):

    """
    Prepare a Camino-style scheme file

    :param subj: utils.Subj
    :param study: utils.Study

    NOTE: this creates a rounded b-value version cuz that's all i need right now. Refactor this function later...

    """

    # we need rounded bvals file. Just call this function. It will check if it needs to create the file
    round_bvals_file(subj,study)


    in_bvals_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvals.rounded")
    if not os.path.isfile(in_bvals_path):
        raise utils.CmdProcError("Missing input:\n"+in_bvals_path)

    in_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvecs.adj")
    if not os.path.isfile(in_bvecs_path):
        in_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvecs")
        if not os.path.isfile(in_bvecs_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvecs_path)


    out_scheme_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.camino.scheme")


    if not os.path.isfile(out_scheme_path):

        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : generate_camino_scheme" )

        cmd = utils.Cmd( "fsl2scheme -bvecfile "+in_bvecs_path+" -bvalfile "+in_bvals_path+" -bscale 1 > "+out_scheme_path )
        cmd.run()







def fsl_bedpostx(subj,study):

    """
    Run FSL bedpostx
    """

    # define working directory
    bpx_work_path = os.path.join(subj.path_append(study, proc_dir),"bpx")


    # check if bedpostx even needs to be run
    if not os.path.isfile(bpx_work_path + ".bedpostX" + os.sep + "dyads1.nii.gz"):


        # define inputs and check
        in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.masked.nii.gz")
        if not os.path.isfile(in_dwi_path):
            in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.masked.nii.gz")
            if not os.path.isfile(in_dwi_path):
                in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.ecc.masked.nii.gz")
                if not os.path.isfile(in_dwi_path):
                    raise utils.CmdProcError("Missing input:\n"+in_dwi_path)

        in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.mask.nii.gz")
        if not os.path.isfile(in_dwi_mask_path):
            in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.mask.nii.gz")
            if not os.path.isfile(in_dwi_mask_path):
                in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.ecc.mask.nii.gz")
                if not os.path.isfile(in_dwi_mask_path):
                    raise utils.CmdProcError("Missing input:\n"+in_dwi_mask_path)

        in_bvals_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvals")
        if not os.path.isfile(in_bvals_path):
            raise utils.CmdProcError("Missing input:\n"+in_bvals_path)

        in_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvecs.adj")
        if not os.path.isfile(in_bvecs_path):
            in_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvecs")
            if not os.path.isfile(in_bvecs_path):
                raise utils.CmdProcError("Missing input:\n"+in_bvecs_path)




        # create working directory
        # try to make path and fail for reasons other than it already exists
        try:
            os.makedirs(bpx_work_path,0o770)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


        # meet FSL file expectations ( in_dwi_path as "data.nii.gz"; brain mask volume in_dwi_mask_path as "nodif_brain_mask.nii.gz", in_bvals_path as 'bvals'; in_bvecs_path as 'bvecs')
        try:

            import shutil as sh
            dict_input_names = {'data.nii.gz':in_dwi_path,'nodif_brain_mask.nii.gz':in_dwi_mask_path,'bvals':in_bvals_path,'bvecs':in_bvecs_path}
            for trg_name, src_path in list(dict_input_names.items()):
                if not os.path.isfile(os.path.join(bpx_work_path,trg_name)):
                    sh.copy(src_path,os.path.join(bpx_work_path,trg_name))
        except:
            raise


        # launch bedpostx
        utils.Cmd('bedpostx '+ bpx_work_path +' -n 2 -w 1 -b 1000 -j 1250 -s 25 -model 2 -c').run()

        # bedpostx command creates bpx for data. This can be deleted manually after bedpostx has run completely but not before. Tricky since we are farming it out to cluster




def extract_topup_b0(subj, study, out_file, b0_indices = ['0', '44']):

    """
    collects single AP PA pair for topup
    because for trio data there aren't dedicated fmaps
    use the first and second run, just hardcoded since it is a one-time application by default but allow other input to accomadate bad acquisitions
    
    """
    in_file_path = subj.path_append(study, proc_dir,  study.label + "." + subj.ursi.short + "." + "dwi.concat.nii.gz")

    # check if already exists
    if not os.path.isfile( out_file ):
        utils.Cmd( "3dcalc -a " + in_file_path + "'[" + ",".join(b0_indices) + "]' -expr 'a' -prefix " + out_file ).run()






def extract_all_b0_to_volume(subj, study, out_file, use_bias_in=False):

    """
    collects all the unweighted images from a raw dwi file an creates a new volume with these only
    NOTE: does not use eddy corrected input here

    :param subj: utils.Subj
    :param study: utils.Study
    :param use_bias_in [True|False]
    :param out_file "Full file path for b0 output volume

    """

    dwi_path          = subj.get_path(study) + "dwi/"
    dwi_file_prefix   = study.label + "." + subj.ursi.short + "." + "dwi.concat"

    if use_bias_in:
        dwi_file_prefix = dwi_file_prefix + ".bc"


    if not os.path.isfile( out_file ):

        if not os.path.isfile(dwi_path + study.label + "." + subj.ursi.short + "." + "dwi.concat.bvecs"):
            raise utils.CmdProcError("Missing bvec file for determining b0 volumes")


        b0_indices = []
        with open(dwi_path + study.label + "." + subj.ursi.short + "." + "dwi.concat.bvecs") as f:
            line_num = 0
            for line in f:
                if line.startswith('0 0 0'):
                    b0_indices.append(str(line_num))
                line_num += 1

        utils.Cmd( "3dcalc -a " + dwi_path + dwi_file_prefix + ".nii.gz'[" + ",".join(b0_indices) + "]' -expr 'a' -prefix " + out_file ).run()





def generate_smt_estimates(subj, study, use_bias_in=False, use_rician_noise_est=True):

    """
    required inputs:
    BS.00638.dwi.concat.bc.edc.eddy_rotated_bvecs (since eddy did motion correction use those bvecs)
    BS.00638.dwi.concat.bc.edc.masked.nii.gz      (to speed loading)
    BS.00638.dwi.concat.bc.edc.mask.nii.gz        (to speed fit)
    BS.00638.dwi.concat.bvals                     (actual bvals)

    default outputs:

    https://github.com/ekaden/smt

    -------------------------
    fitmicrodt
    Microscopic diffusion tensor
    This microscopic model estimates the microscopic fractional anisotropy and other microscopic diffusion coefficients in tissue with complex directional architecture as in the brain.
    -------------------------
    Longitudinal microscopic diffusivity
    Transverse microscopic diffusivity
    Microscopic fractional anisotropy
    Microscopic mean diffusivity
    Zero b-value image


    -------------------------
    fitmcmicro
    Multi-compartment microscopic diffusion
    This model provides estimates of microscopic features specific to the intra- and extra-neurite compartments unconfounded by the effects of the potentially complex fibre orientation distribution.
    -------------------------
    Intra-neurite volume fraction
    Intrinsic diffusivity
    Extra-neurite transverse microscopic diffusivity
    Extra-neurite microscopic mean diffusivity
    Zero b-value image

    """


    # define input directory
    dwi_path = subj.path_append(study, proc_dir)
    # define output directory
    smt_path = os.path.join(subj.path_append(study, proc_dir),"smt")



    in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.masked.nii.gz")
    in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.edc.mask.nii.gz")
    in_bvecs_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.eddy_rotated_bvecs")
    in_bvals_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bvals")
    out_file_prefix = study.label + "." + subj.ursi.short + "." + "dwi.concat.edc.masked"

    if use_bias_in:
        in_dwi_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.masked.nii.gz")
        in_dwi_mask_path = subj.path_append(study, proc_dir) + subj.id_fn(study,"dwi.concat.bc.edc.mask.nii.gz")
        out_file_prefix = study.label + "." + subj.ursi.short + "." + "dwi.concat.bc.edc.masked"




    if not os.path.isfile(in_dwi_path):
        raise utils.CmdProcError("Missing input:\n"+in_dwi_path)

    if not os.path.isfile(in_dwi_mask_path):
        raise utils.CmdProcError("Missing input:\n"+in_dwi_mask_path)

    if not os.path.isfile(in_bvals_path):
        raise utils.CmdProcError("Missing input:\n"+in_bvals_path)

    if not os.path.isfile(in_bvecs_path):
        raise utils.CmdProcError("Missing input:\n"+in_bvecs_path)



    # create smt directory
    try:
        os.makedirs(smt_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


    # estimate rician noise
    if use_rician_noise_est:

        # desired output
        o_file_rician_noise_est = os.path.join(smt_path, "rician_noise_est.nii.gz")
        if not os.path.isfile(o_file_rician_noise_est):

            # create a b0 only volume for noise estimation
            # using edc input in this case
            b0_filepath = os.path.join(smt_path, out_file_prefix + ".b0.nii.gz")
            if not os.path.isfile( b0_filepath ):

                b0_indices = []
                with open(dwi_path + study.label + "." + subj.ursi.short + "." + "dwi.concat.bvecs") as f:
                    line_num = 0
                    for line in f:
                        if line.startswith('0 0 0'):
                            b0_indices.append(str(line_num))
                        line_num += 1

                utils.Cmd( "3dcalc -a " + in_dwi_path + "'[" + ",".join(b0_indices) + "]' -expr 'a' -prefix " + b0_filepath ).run()


            tmp_snr_file = os.path.join(smt_path, "rician_snr_est.nii")

            if not os.path.isfile( tmp_snr_file ):
                utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : generating rician noise estimate" )
                utils.Cmd("ricianfit --mask "+in_dwi_mask_path+" "+b0_filepath+" " +tmp_snr_file ).run()
                #utils.Cmd("ricianfit "+b0_filepath+" " +os.path.join(smt_path, "rician_snr_est") ).run()

            # two volume file: 0=signal; 1=noise
            if os.path.isfile( tmp_snr_file ):
                # split out the noise volume from the output
                # Don't use AFNI as it will "correct" float errors and create zeros in the noise estimate volume which will crash fitmcmicro cuz it doesn't want to divide by zero
                utils.Cmd("fslsplit "+tmp_snr_file+" "+os.path.join(smt_path, "rician_snr")+" -t ").run()
                os.rename( os.path.join(smt_path, "rician_snr0000.nii.gz"), os.path.join(smt_path, "rician_signal_est.nii.gz") )
                os.rename( os.path.join(smt_path, "rician_snr0001.nii.gz"), o_file_rician_noise_est)

                # clean up
                utils.purge_paths(b0_filepath, tmp_snr_file)



    # fit Multi-compartment microscopic diffusion
    #Gaussian
    final_o_labels_gaus = ['IntraNeuriteVolFrac-Gaus','IntrinsicDiff-Gaus','ExtraNeuriteTransDiff-Gaus','ExtraNeuriteMeanDiff-Gaus','mcmicro-ZeroB-Gaus']
    final_o_file = os.path.join(smt_path, final_o_labels_gaus[-1]+".nii.gz")

    if not os.path.isfile(final_o_file):

        tmp_fit_file = os.path.join(smt_path, "mcmicro_est_gaussian.nii")
        if not os.path.isfile(tmp_fit_file):
            utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Fit Multi-compartment microscopic diffusion : Gaussian" )
            utils.Cmd("fitmcmicro --bvals "+in_bvals_path+" --bvecs "+in_bvecs_path+" --mask "+in_dwi_mask_path+" "+in_dwi_path+" "+tmp_fit_file ).run()


        if os.path.isfile(tmp_fit_file):
            #Intra-neurite volume fraction
            #Intrinsic diffusivity
            #Extra-neurite transverse microscopic diffusivity
            #Extra-neurite microscopic mean diffusivity
            #Zero b-value image
            for brik_indx,brik_label in enumerate(final_o_labels_gaus):
                o_file = os.path.join(smt_path, brik_label+".nii.gz")
                if not os.path.isfile(o_file):
                    utils.Cmd("3dcalc -a "+tmp_fit_file+"'["+str(brik_indx)+"]' -expr 'a' -prefix "+o_file).run()

        utils.purge_paths(tmp_fit_file)





    # fit Multi-compartment microscopic diffusion
    #Rician
    if use_rician_noise_est:

        final_o_labels_ric = ['IntraNeuriteVolFrac-Ric','IntrinsicDiff-Ric','ExtraNeuriteTransDiff-Ric','ExtraNeuriteMeanDiff-Ric','mcmicro-ZeroB-Ric']
        final_o_file = os.path.join(smt_path, final_o_labels_ric[-1]+".nii.gz")

        if not os.path.isfile(final_o_file):

            tmp_fit_file = os.path.join(smt_path, "mcmicro_est_rician.nii")
            if not os.path.isfile(tmp_fit_file):
                utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : DWI : Fit Multi-compartment microscopic diffusion : Rician" )
                utils.Cmd("fitmcmicro --bvals "+in_bvals_path+" --bvecs "+in_bvecs_path+" --rician "+o_file_rician_noise_est+" --mask "+in_dwi_mask_path+" "+in_dwi_path+" "+tmp_fit_file ).run()


            if os.path.isfile(tmp_fit_file):
                #Intra-neurite volume fraction
                #Intrinsic diffusivity
                #Extra-neurite transverse microscopic diffusivity
                #Extra-neurite microscopic mean diffusivity
                #Zero b-value image
                for brik_indx,brik_label in enumerate(final_o_labels_ric):
                    o_file = os.path.join(smt_path, brik_label+".nii.gz")
                    if not os.path.isfile(o_file):
                        utils.Cmd("3dcalc -a "+tmp_fit_file+"'["+str(brik_indx)+"]' -expr 'a' -prefix "+o_file).run()


            utils.purge_paths(tmp_fit_file)


    ## label subriks in output or  spllit output to separate volumes

    ## fit Microscopic diffusion tensor
    #Gaussian
    #fitmicrodt --bvals BS.00638.dwi.concat.bvals.rounded --bvecs BS.00638.dwi.concat.bc.edc.eddy_rotated_bvecs --mask BS.00638.dwi.concat.bc.edc.mask.nii.gz BS.00638.dwi.concat.bc.edc.masked.nii.gz fitmicrodt_gaussian_estimates
    #Rician
    #xxxx fitmicrodt --bvals BS.00638.dwi.concat.bvals.rounded --bvecs BS.00638.dwi.concat.bc.edc.eddy_rotated_bvecs --mask BS.00638.dwi.concat.bc.edc.mask.nii.gz BS.00638.dwi.concat.bc.edc.masked.nii.gz fitmicrodt_gaussian_estimates

