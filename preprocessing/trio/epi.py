"""@package epi
Preprocessing module for functional imaging

"""

from mayerlab import utils    #
import os                #
import errno            #    for path check and create





def epi2tlrc_affine(subject,study,run_type,cost_func="lpc",move="",opts="",purge=False):
    """
    affine alignment of epi to T1 and TLRC space

        #!/usr/bin/env python
    change to
        #!/export/research/analysis/human/amayer/shared/apps/python/bin/python
    in align_epi_anat.py
    to avoid import error
    """

    if purge or not os.path.isfile(subject.get_transform_path(study) + study.label + "." + subject.ursi.short + "." + run_type + "_epi2tlrc.aff12.1D"):

        #    for directory clean up
        import shutil


        workspace_path = subject.get_path(study) + "anat/regTemp." + run_type

        #    create alignment workspace
        try:
            os.makedirs(workspace_path,0o770)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        #    go in there cuz afni script is too dumb to use full paths correctly
        current_path = os.getcwd()
        os.chdir( workspace_path )


        ##    copying over working files
        cmd = utils.Cmd("cp " + subject.get_path(study) + "anat/" + study.label + "." + subject.ursi.short + ".anat.sksp+* " + workspace_path + "/" )
        cmd.run()

        cmd = utils.Cmd("cp " + subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".BASEIMG_sksp+orig.* " + workspace_path + "/" )
        cmd.run()


        # this parameter comes to use without a hyphen to keep it from being interpreted as a command. Drop hyphen on this string if it user passed in an argument
        if len(move):
               move = "-"+move

        #    do alignment, sadly afni script not made callable from python
        cmd = utils.Cmd("align_epi_anat.py -epi2anat -anat " + workspace_path + "/" + study.label + "." + subject.ursi.short + ".anat.sksp+orig -epi " + workspace_path + "/" + study.label + "." + subject.ursi.short + "." + run_type + ".BASEIMG_sksp+orig -epi_base 0 -volreg off -epi_strip None -deoblique on -anat_has_skull no -tshift off -align_centers yes -giant_move -cost " + cost_func + " " + move + " -tlrc_apar " + workspace_path + "/" + study.label + "." + subject.ursi.short + ".anat.sksp+tlrc -Allineate_opts " + opts + " -source_automask+4 ")
        cmd.run()


        #    rename matrices
        cmd = utils.Cmd("cat_matvec -ONELINE " + workspace_path + "/" + study.label + "." + subject.ursi.short + "." + run_type + ".BASEIMG_sksp_al+orig::ALLINEATE_MATVEC_B2S_000000 > " + subject.get_transform_path(study) + study.label + "." + subject.ursi.short + "." + run_type + "_epi2T1.aff12.1D" )
        cmd.run()

        cmd = utils.Cmd("cp " + workspace_path + "/" + study.label + "." + subject.ursi.short + "." + run_type + ".BASEIMG_sksp_al_tlrc_mat.aff12.1D " + subject.get_transform_path(study) + study.label + "." + subject.ursi.short + "." + run_type+ "_epi2tlrc.aff12.1D" )
        cmd.run()


        #    clean up
        os.chdir( current_path )
        shutil.rmtree(workspace_path, ignore_errors=True)

        #    since affine align has been updated we need to purge final run alignment if it exists
        utils.purge_paths(subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".preproc.3x3x3+tlrc.*")







def qa_epi2tlrc_affine(subject,study,run_type):

    """
    check final affine/qwarp alignment of epi to T1 and TLRC space
    first check affine, then check non-linear results
    """

    #    try doing this in machine tmp since this is all temporary output
    os.chdir("/tmp")


    utils.prettyOut(run_type + " : EPI affine alignment check")

    if not os.path.isfile("/tmp/VERIFY." + run_type + ".tmp+tlrc.HEAD"):

        cmd = utils.Cmd("3dAllineate " +
                " -master " + subject.get_path(study) + "anat/" + study.label + "." + subject.ursi.short + ".anat.sksp+tlrc " +
                " -1Dmatrix_apply " + subject.get_transform_path(study) + study.label + "." + subject.ursi.short + "." + run_type + "_epi2tlrc.aff12.1D" +
                " -input " + subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".BASEIMG_sksp+orig " +
                " -floatize" +
                " -interp trilinear" +
                " -final wsinc5"+
                " -prefix /tmp/VERIFY." + run_type + ".tmp" )

        cmd.run()


    ##    copying over working files
    cmd = utils.Cmd("cp " + subject.get_path(study) + "anat/" + study.label + "." + subject.ursi.short + ".anat.sksp+tlrc* " + "/tmp/" )
    cmd.run()


    cmd = utils.Cmd("3dedge3 -overwrite -prefix /tmp/VERIFY." + run_type + ".tmp.edge -input  /tmp/VERIFY." + run_type + ".tmp+tlrc.HEAD" )
    cmd.run()


    # scale the volume for better contrast
    cmd = utils.Cmd("@ScaleVolume -input /tmp/VERIFY." + run_type + ".tmp.edge+tlrc -prefix /tmp/VERIFY." + run_type + ".tmp.edge.sc -perc_clip 2 98 -val_clip 6 249")
    cmd.run()

    #    Don't let afni session detach or it will close before user can see images
    os.environ["AFNI_DETACH"] = "NO"

    cmd = utils.Cmd('afni -com "SWITCH_UNDERLAY ' + study.label + '.' + subject.ursi.short + '.anat.sksp" -com "SWITCH_OVERLAY VERIFY.' + run_type + '.tmp.edge.sc" -com "SET_XHAIRS OFF" -com "SEE_OVERLAY +" -com "SET_PBAR_SIGN A.+" -com "SET_FUNC_AUTORANGE A.+" -com "SET_THRESHOLD A.40 2" -com "SET_SPM_XYZ 6 12 14" ')
    cmd.run()

    #    return this setting to lab default
    os.environ["AFNI_DETACH"] = "YES"


    #    clean up
    utils.purge_paths("/tmp/" + study.label + "." +  subject.ursi.short + ".anat.*","/tmp/VERIFY." + run_type + ".*")







def epi2tlrc_nonlinear(subject,study,run_type,regrid=3,base_suffix="dist_corr"):

    """
    apply linear and non-linear warp to volume to prepare for level1

    """


    if not os.path.isfile(subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".preproc.3x3x3+tlrc.HEAD"):

        #    make sure required input exists
        if not os.path.isfile(subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".preproc."+base_suffix+"+orig.HEAD"):
            raise utils.CmdProcError("Missing concatenated run volume input: " + subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".preproc."+base_suffix+"+orig.HEAD +++ Make sure EPI preprocessing has completed without error.")



        utils.prettyOut(run_type + " : Applying non-linear alignment to preproc volume for level1. Patience. This could take awhile...")

        src_path            = subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".preproc."+base_suffix+"+orig"
        nl_warp_src        = subject.get_transform_path(study) + study.label+"."+subject.ursi.short+".anat.sksp.qwarp_WARP1+tlrc"
        t1_to_tlrc_mat    = subject.get_transform_path(study) + study.label+"."+subject.ursi.short+".T1_2tlrc.aff3x4.1D"
        epi_2_T1_aff_mat    = subject.get_transform_path(study) + study.label+"."+subject.ursi.short+"."+run_type+"_epi2T1.aff12.1D"


        cmd = utils.Cmd("3dNwarpApply " +
                " -source " + src_path  +
                " -nwarp '" + nl_warp_src + " " + t1_to_tlrc_mat + " " + epi_2_T1_aff_mat + "'" +
                " -master " + subject.get_path(study) + "anat/" + study.label+"."+subject.ursi.short+".anat.sksp.qwarp+tlrc" +
                " -interp wsinc5" +
                " -newgrid " + str(regrid) +
                " -prefix " + subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".preproc.3x3x3" )

        cmd.run()






def qa_epi2tlrc_nonlinear(subject,study,run_type,base_suffix="dist_corr"):

    """
    check final affine/qwarp alignment of epi to T1 and TLRC space
    first check affine, then check non-linear results
    """

    #    try doing this in machine tmp since this is all temporary output
    os.chdir("/tmp")




    ##    copying over working files from anat
    cmd = utils.Cmd("cp " + subject.get_path(study) + "anat/" + study.label+"."+subject.ursi.short+".anat.sksp.qwarp+tlrc*" + " /tmp/" )
    cmd.run()


    #    define volumes to check ( a label for screen and input path
    #ref_vol_dict = { '2 epi data sample volume':subject.get_path(study) + "REG_VOLUMES/"  + study.label+"."+subject.ursi.short+"."+run_type+".preproc."+base_suffix+"+orig[0]" , '1 reference volume':subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".BASEIMG_sksp+orig" }
    ref_vol_dict = { '1 reference volume':subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type + ".BASEIMG_sksp+orig" , '2 epi data sample volume':subject.get_path(study) + "REG_VOLUMES/"  + study.label+"."+subject.ursi.short+"."+run_type+".preproc."+base_suffix+"+orig[0]" }
    #ref_vol_dict = { 'epi data sample volume':subject.get_path(study) + "REG_VOLUMES/"  + study.label+"."+subject.ursi.short+"."+run_type+".preproc."+base_suffix+"+orig[0]" }
    
    # loop over volumes to check
    for src_label, src_path in list(ref_vol_dict.items()):

        #    message to user
        utils.prettyOut(run_type + " : " + src_label + " : EPI non-linear alignment check")
        utils.prettyOut("Task (aligned to)-> SBREF (aligned to)-> T1 (aligned to)-> template")

        #    create verify volume if not there
        if not os.path.isfile("/tmp/" + study.label + "." + subject.ursi.short + "." + run_type + ".VERIFY_sksp+tlrc.HEAD"):

            #    do SBRef due to nice contrast for visual check
            nl_warp_src        = subject.get_transform_path(study) + study.label+"."+subject.ursi.short+".anat.sksp.qwarp_WARP1+tlrc"
            t1_to_tlrc_mat    = subject.get_transform_path(study) + study.label+"."+subject.ursi.short+".T1_2tlrc.aff3x4.1D"
            epi_2_T1_aff_mat    = subject.get_transform_path(study) + study.label+"."+subject.ursi.short+"."+run_type+"_epi2T1.aff12.1D"

            #    the reference volume output will be aligned to
            #    anat.sksp.qwarp+tlrc
            #    which is non-linear warp to TT_N27
            #    so reference volume output should look good over that
            cmd = utils.Cmd("3dNwarpApply " +
                    " -source " + src_path  +
                    " -nwarp '" + nl_warp_src + " " + t1_to_tlrc_mat + " " + epi_2_T1_aff_mat + "'" +
                    " -master " + subject.get_path(study) + "anat/" + study.label+"."+subject.ursi.short+".anat.sksp.qwarp+tlrc" +
                    " -interp wsinc5" +
                    " -prefix " + "/tmp/" + study.label + "." + subject.ursi.short + "." + run_type + ".VERIFY_sksp+tlrc.HEAD" )

            cmd.run()



        #    edge for ease of viewing
        cmd = utils.Cmd("3dedge3 -overwrite -prefix /tmp/" + study.label + "." + subject.ursi.short + "." + run_type + ".VERIFY_sksp.edge -input  /tmp/" + study.label + "." + subject.ursi.short + "." + run_type + ".VERIFY_sksp+tlrc.HEAD" )
        cmd.run()


        # scale the volume for better contrast
        cmd = utils.Cmd("@ScaleVolume -input /tmp/" + study.label + "." + subject.ursi.short + "." + run_type + ".VERIFY_sksp.edge+tlrc -prefix /tmp/" + study.label + "." + subject.ursi.short + "." + run_type + ".VERIFY_sksp.sc -perc_clip 2 98 -val_clip 6 249")
        cmd.run()


        #    Don't let afni session detach or it will close before user can see images
        os.environ["AFNI_DETACH"] = "NO"

        cmd = utils.Cmd('afni -com "SWITCH_UNDERLAY ' + study.label + '.' + subject.ursi.short + '.anat.sksp.qwarp" -com "SWITCH_OVERLAY ' + study.label + '.' + subject.ursi.short + '.' + run_type + '.VERIFY_sksp.sc" -com "SET_XHAIRS OFF" -com "SEE_OVERLAY +" -com "SET_PBAR_SIGN A.+" -com "SET_FUNC_AUTORANGE A.+" -com "SET_THRESHOLD A.40 2" -com "SET_SPM_XYZ 9 21 18" ')
        cmd.run()

        #    return this setting to lab default
        os.environ["AFNI_DETACH"] = "YES"

        #    clean up verify volume
        utils.purge_paths("/tmp/" + study.label + "." + subject.ursi.short + "." + run_type + ".VERIFY_sksp*")



    #    clean up anat at the end
    utils.purge_paths("/tmp/" + study.label + "." +  subject.ursi.short + ".anat.sksp.qwarp*")





def to3d(subj,study,run_label):

    """
    create volume from dicom

    """
    run_settings = study.get_run_label_settings(run_label)
    out_path = subj.get_path(study) + run_label + "/"
    out_file = study.label + "." + subj.ursi.short + "." + run_label

    #    try to make output path and fail for reasons other than it already exists
    try:
        os.makedirs(out_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    
    if not os.path.isfile(out_path + out_file + "+orig.HEAD"):
        
        import glob

        utils.prettyOut(run_label + " : Create volume")

        #    need dicom path from _run_map_anat_
        dcm_path = subj.get_dcm_path(study,"epi",run_label)
        if not dcm_path:
            raise utils.CmdProcError("Dicom path not found")
            
        #    let's do a file count just so we can provide a more meaningful error than what to3d will say.
        dcm_iter = glob.iglob(dcm_path + "/*.dcm")
        dcm_iter_len = len(list(dcm_iter))

        if(int(dcm_iter_len) != int(run_settings['time_images']) ):
            raise utils.CmdProcError("Not enough dicom files at {}\n Expected {} but found {}".format(dcm_path + "/*.dcm",run_settings['time_images'],dcm_iter_len) )
        
        
        # to get correct slice timing and check dicom files use Dimon to get file list
        dcm_org_file_list = os.path.join(out_path, out_file+'.dicom_file_list')
        if not os.path.isfile(dcm_org_file_list):
            # run Dimon
            utils.Cmd("Dimon -infile_pattern '" + dcm_path + "/*.dcm' -dicom_org -save_file_list "+dcm_org_file_list+" -quit").run()            
            
        utils.Cmd("to3d -time:zt " + run_settings['num_slices'] + " " + run_settings['time_images'] + " " + run_settings['TR'] + " FROM_IMAGE -save_outliers " + out_path + out_file + "_outliers -epan -session " + out_path + " -orient RAI -prefix " + out_file + " -@ < " + dcm_org_file_list).run()






def to3d_SBRef(subj, study, run_type):

    """
    create volume from dicom

    """

    #    we'll put output into RUN1 directory
    run_label    = run_type + "1"
    out_path        = subj.get_path(study) + run_label + "/"
    #    filename
    out_file        = study.label + "." + subj.ursi.short + "." + run_label + ".SBRef"

    #    try to make output path and fail for reasons other than it already exists
    try:
        os.makedirs(out_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    if not os.path.isfile(out_path + out_file + "+orig.HEAD"):

        #    need dicom path from _run_map_anat_
        try:
            dcm_path = subj.get_dcm_path(study,"epi",run_label)
        except Exception as exception:
            raise exception
            
        if not dcm_path:
            raise utils.CmdProcError("Could not find dicom path for key: " + run_label + ". Check that run_map file contains a line key: "+ run_label)
            
        #    pull elements from path so that i derive SBREF dicom src from run type
        dcm_path_elements = utils.parse_dicom_path(dcm_path)
        if not len(dcm_path_elements['series']):
            raise utils.CmdProcError("Could not parse found dicom path for " + run_label)
            
        #    because SBREF could be rerun we have to either take the last one or be more specific
        #    i choose to be precise. the path below takes the run series from the run_map listing and decrements it
        #    to specify the series number that should be mated to this run by virtue of how the scanner sequence is programmed
        #    here the length of series is fed into the format specifiers as replacement field as a small defense against future change
        try:
            sbref_path = dcm_path_elements['directory'] + "/" + dcm_path_elements['sequence'] + "_SBRef_" + "{:0>{}}".format( int(dcm_path_elements['series'])-1, len(dcm_path_elements['series']) ) + "/*.dcm"
        except ValueError as exception:
            raise exception

        utils.prettyOut(run_type + " : Create SBRef volume")
        utils.Cmd("to3d -epan -session " + out_path + " -orient RAI -prefix " + out_file + " " + sbref_path).run()



def distcorr_to_SBRef(subj,study,run_type):

    """
    Rigid alignment of distortion correction spin echo to SBREF gradient echo
    This is so that distortion correction estimate and tissue match.
    """

    #    # ls rigid alignment of original dist correction volume to SBREF for transform matrix
    #    # RADIOLOGICAL
    #    flirt -in distortion_corr_32ch_APPA.nii.gz -ref LAPTOP.70647.CVR1.SBRef.nii.gz -omat flirt_distcorr_to_SBRef.mat -bins 256 -cost leastsq -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6 -interp trilinear
    #
    #    # create temp version of field coefficients
    #    cp topup_results_APPA_fieldcoef.nii.gz NEWHD_topup_results_APPA_fieldcoef.nii.gz
    #
    #    # alter header of field coefficients temp to match field map
    #    fslcpgeom topup_field_APPA.nii.gz NEWHD_topup_results_APPA_fieldcoef.nii.gz
    #
    #    # now apply transform to that version
    #    flirt -in NEWHD_topup_results_APPA_fieldcoef.nii.gz -applyxfm -init flirt_distcorr_to_SBRef.mat -out topup_results_APPA_fieldcoef_SBREF -paddingsize 0.0 -interp trilinear -ref NEWHD_topup_results_APPA_fieldcoef.nii.gz
    #
    #    # reinstate original header of transformed output.
    #    # NOTE: field map output from topup is RADIOLOGICAL
    #    #       field coefficients from topup is NEUROLOGICAL
    #    fslcpgeom topup_results_APPA_fieldcoef.nii.gz topup_results_APPA_fieldcoef_SBREF.nii.gz
    #
    #    rm NEWHD_topup_results_APPA_fieldcoef.nii.gz


    run_label        = run_type + "1"
    run_in_path      = subj.get_path(study) + run_label + "/"
    dist_in_path     = os.path.join(subj.get_transform_path(study), "distortion_corr_32ch")

    ref_nii_file     = os.path.join(run_in_path, study.label + "." + subj.ursi.short + "." + run_label + ".SBRef.nii.gz")
    ref_afni_file    = os.path.join(run_in_path, study.label + "." + subj.ursi.short + "." + run_label + ".SBRef+orig")

    dc_i_file        = os.path.join(dist_in_path, "distortion_corr_32ch_APPA.nii.gz")
    dc_i_mov_file    = os.path.join(dist_in_path, "topup_results_APPA_movpar.txt")
    dc_i_field_file  = os.path.join(dist_in_path, "topup_field_APPA.nii.gz")
    dc_i_coef_file   = os.path.join(dist_in_path, "topup_results_APPA_fieldcoef.nii.gz")
    dc_i_coef_tmp    = os.path.join(dist_in_path, "NEWHD_topup_results_APPA_fieldcoef.nii.gz")

    # NOTE: applytopup doesn't like periods in filename
    dc_o_coef_file   = os.path.join(run_in_path, "topup_results_APPA_SBRef_fieldcoef.nii.gz")
    dc_o_mov_file    = os.path.join(run_in_path, "topup_results_APPA_SBRef_movpar.txt")
    mat_o_file       = os.path.join(subj.get_transform_path(study), study.label + "." + subj.ursi.short + "." + run_label + ".flirt_distcorr_to_SBRef_" + run_label + ".mat")


    # do we already have what we need?
    if not os.path.isfile(dc_o_coef_file):

        utils.prettyOut("Align distortion correction spin echo to SBREF gradient echo " + ref_afni_file)

        # verify inputs
        if not os.path.isfile( ref_afni_file + ".HEAD"):
            raise utils.CmdProcError("Missing " + ref_afni_file + ".HEAD")
        if not os.path.isfile( dc_i_file ):
            raise utils.CmdProcError("Missing " + dc_i_file)

        # NIFTI format for FSL
        if not os.path.isfile(ref_nii_file):
            utils.Cmd("3dresample -orient RPI -inset " + ref_afni_file + " -prefix " + ref_nii_file).run()

        # Get transformation matrix. rigid alignment of original dist correction volume to SBREF for transform matrix
        if not os.path.isfile(mat_o_file):
            utils.Cmd("flirt -in "+dc_i_file+" -ref "+ref_nii_file+" -omat "+mat_o_file+" -bins 256 -cost leastsq -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6 -interp trilinear").run()

        # create temp version of field coefficients
        import shutil
        shutil.copy(dc_i_coef_file, dc_i_coef_tmp)

        # alter header of field coefficients temp to match field map
        utils.Cmd("fslcpgeom "+dc_i_field_file+" "+dc_i_coef_tmp).run()

        # now apply transform to that version
        utils.Cmd("flirt -in "+dc_i_coef_tmp+" -applyxfm -init "+mat_o_file+" -out "+dc_o_coef_file+" -paddingsize 0.0 -interp trilinear -ref " + dc_i_coef_tmp).run()

        # reinstate original header of transformed output.
        # NOTE: field map output from topup is RADIOLOGICAL
        #       field coefficients from topup is NEUROLOGICAL
        utils.Cmd("fslcpgeom "+dc_i_coef_file+" "+dc_o_coef_file).run()

        # applytopup expects mov file at same location as coef file
        shutil.copy(dc_i_mov_file, dc_o_mov_file)

        # remove temp
        utils.purge_paths(dc_i_coef_tmp)






def despike(subj,study,run_label,ignore=1):

    """
    despike

    """

    out_path = subj.get_path(study) + run_label + "/"
    out_file = out_path + study.label + "." + subj.ursi.short + "." + run_label + ".dspk"

    #    make sure required input exists
    if not os.path.isfile( out_path + study.label + "." + subj.ursi.short + "." + run_label + "+orig.HEAD" ):
        raise utils.CmdProcError("Missing raw volume input")

    #    only run if output does not already exist
    if not os.path.isfile( out_file + "+orig.HEAD"):

        utils.prettyOut(run_label + " : Despike")

        cmd = utils.Cmd("3dDespike -NEW -ignore " + str(ignore) + " -nomask -cut 4 4.5 -prefix "+out_file+" " + out_path + study.label + "." + subj.ursi.short + "." + run_label + "+orig" )
        cmd.run()





def timeshift(subj,study,run_label,ignore=1):

    """
    timeshift

    """

    out_path = subj.get_path(study) + run_label + "/"
    out_file = out_path + study.label + "." + subj.ursi.short + "." + run_label + ".tshft"

    #    only run if output does not already exist
    if not os.path.isfile(out_file + "+orig.HEAD"):

        utils.prettyOut(run_label + " : Timeshift")

        cmd = utils.Cmd("3dTshift -ignore "+str(ignore)+ " -tzero 0 -prefix "+out_file+" " + out_path + study.label + "." + subj.ursi.short + "." + run_label + ".dspk+orig" )
        cmd.run()






def imReg_2d(subj,study,run_label,base_suffix=".SBRef"):

    """
    2d registration to SBREF base file

    """

    out_path     = subj.get_path(study) + run_label + "/"
    out_file     = out_path + study.label + "." + subj.ursi.short + "." + run_label + ".2dreg"    # used twice; 1st as volume prefix 2nd as transform prefix
    
    in_file      = out_path + study.label + "." + subj.ursi.short + "." + run_label + ".tshft+orig"
    if not os.path.isfile(in_file + "+orig.HEAD"):
        in_file      = out_path + study.label + "." + subj.ursi.short + "." + run_label + ".dspk+orig"
        
    run_type        = run_label.rstrip("0123456789")
    #    base file should be from run1 of this run type but you can specify which of those in call
    base_file     = subj.get_path(study) + run_type + "1/" + study.label + "." + subj.ursi.short + "." + run_type + "1" + base_suffix + "+orig"


    #    only run if output does not already exist
    if not os.path.isfile(out_file + "+orig.HEAD"):

        utils.prettyOut( run_label + " : DATA IS REGISTERED IN 2D ")

        cmd = utils.Cmd("2dImReg -input " + in_file + " -basefile " + base_file + " -base 0 -prefix " + out_file + " -dprefix " + out_file )
        cmd.run()






def imReg_3d(subj,study,run_label,base_suffix=".SBRef"):

    """
    3d registration of processed data to SBREF base file to correct motion

    This outputs parameters used as motion regressors. See create_run_motion_regressors()
    ".Movement.Regressor.mc.1D"):
    """

    out_path     = subj.get_path(study) + run_label + "/"
    out_file     = out_path + study.label + "." + subj.ursi.short + "." + run_label + ".3dreg"
    in_file      = out_path + study.label + "." + subj.ursi.short + "." + run_label + ".2dreg+orig"
    run_type        = run_label.rstrip("0123456789")
    base_file     = subj.get_path(study) + run_type + "1/" + study.label + "." + subj.ursi.short + "." + run_type + "1" + base_suffix + "+orig"

    #    only run if output does not already exist
    if not os.path.isfile(out_file + "+orig.HEAD"):

        utils.prettyOut( run_label + " : DATA IS REGISTERED IN 3D ")
        cmd = utils.Cmd("3dvolreg -dfile " + out_file + "log.mc -1Dmatrix_save " + out_path + study.label + "." + subj.ursi.short + "." + run_label + ".mat.mc -prefix " + out_file + " -base " + base_file + " " + in_file )
        cmd.run()




        



def estimate_run_motion(subj,study,run_label,base_suffix=".SBRef"):

    """
    Estimate movement parameters
    Register the raw, unprocessed run data to SBRef to estimate motion estimates.

    This is used for motion quantification. See create_run_motion_regressors()
    ".Movement.Regressor.1D"):

    """

    out_path     = subj.get_path(study) + run_label + "/"
    out_file     = out_path + study.label + "." + subj.ursi.short + "." + run_label + ".3dreglog.raw"
    run_type        = run_label.rstrip("0123456789")
    base_file     = subj.get_path(study) + run_type + "1/" + study.label + "." + subj.ursi.short + "." + run_type + "1" + base_suffix + "+orig"

    #    only run if output does not already exist
    if not os.path.isfile(out_file):

        utils.prettyOut(run_label + " : Estimate run motion parameters ")
        cmd = utils.Cmd("3dvolreg -dfile " + out_file + " -prefix NULL -base 0 -base " + base_file + " " + out_path + study.label + "." + subj.ursi.short + "." + run_label + "+orig" )
        cmd.run()





def create_motion_deriv_file(subj, study, run_label):

    """
    Create derivative of motion regressor file

    """
    
    #    run_label indicates where processing data is
    #    run_type indicates final processing output
    run_type = run_label.rstrip("0123456789")
    run_type_settings = study.get_run_type_settings(run_type)

    
    i_file = os.path.join( study.movement_path, subj.visit, subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.1D" )    
    o_file = os.path.join( study.movement_path, subj.visit, subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D" )

    #    generate motion derivative regessors if they don't exist
    if not os.path.isfile( o_file ):

        utils.Cmd("1d_tool.py"
                " -infile " + i_file +
                " -set_nruns " + run_type_settings['run_count'] + " "
                " -derivative"
                " -write " + o_file ).run()






def concate_runs(subj,study,run_type,base=1):

    """
    glue together preprocessed runs eliminating anything before base image

    """

    import glob

    out_path = subj.get_path(study) + "REG_VOLUMES/"
    out_file = out_path + study.label +"." + subj.ursi.short + "." + run_type + ".preproc.uncorrected"
    input_vols = glob.glob(subj.get_path(study) + run_type + "?/" + study.label + "." + subj.ursi.short + ".*" + ".3dreg+orig.HEAD")


    #    try to make output path and fail for reasons other than it already exists
    try:
        os.makedirs(out_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise



    if not os.path.isfile(out_file + "+orig.HEAD"):

        utils.prettyOut(run_type + " Concatenate " + str(len(input_vols)) + " runs")

        cmd = utils.Cmd("3dTcat -prefix " + out_file + " " + "'[1-$]' ".join(sorted(input_vols)) + "'[1-$]'" )
        cmd.run()







def generate_run_mask(subj,study,run_type,base_suffix=".SBRef.dist_corr"):

    """
    create mask on base image using final concatenated data for this run type

    """

    out_path    = subj.get_path(study) + "REG_VOLUMES/"

    in_vol     = subj.get_path(study) + run_type.upper() + "1/" + study.label + "." + subj.ursi.short + "." + run_type.upper() + "1" + base_suffix


    if not os.path.isfile(in_vol + "+orig.HEAD"):
        raise utils.CmdProcError("Missing input for masking")

    if not os.path.isfile(out_path + study.label + "." + subj.ursi.short + "." + run_type.upper() + ".BASEIMG_sksp+orig.HEAD"):

        utils.prettyOut(run_type + " : Base run mask")

        cmd = utils.Cmd("3dAutomask -apply_prefix " + out_path + study.label + "." + subj.ursi.short + "." + run_type.upper() + ".BASEIMG_sksp " + in_vol + "+orig'[0]'" )
        cmd.run()








def estimate_field_distortion(subject,study):

    """
    Uses flipped phase encoding acquisitions to estimate b0 field distortions for this subject for this session.

    """

    import glob

    # look for distortion acqusitions in the session in which the epi was acquired
    # will fail if not EPI runs acquired
    epi_run_list = subject.get_run_list(study,mri_type="epi",run_type=None,mapped=True)
    if not len(epi_run_list):
        raise Exception("No EPI runs")

    # session path taken from 1st mapped epi run
    session_path = subject.get_session_path(study,"epi",epi_run_list[0])   
    # look for distortion acquisition at the path
    path_found_list = glob.glob(session_path + "/distortion_corr_32ch_ap_0*")
    
    if not len(path_found_list):
        # look instead at the mprage session. Though this makes little sense if there are no EPI runs
        session_path = subject.get_session_path(study,"anat",'MPRAGE1')
        path_found_list = glob.glob(session_path + "/distortion_corr_32ch_ap_0*")

        if not len(path_found_list):
            raise Exception("No session path for distortion_corr_32ch")
    
    
    
    
    #print( glob.glob(epi_session_path + "/distortion_corr_32ch_ap_0*") )
    # get last run of
    path_ap = glob.glob(session_path + "/distortion_corr_32ch_ap_0*")[-1]
    # get last run of 
    path_pa = glob.glob(session_path + "/distortion_corr_32ch_pa_0*")[-1]

    #    try to make output path and fail for reasons other than it already exists
    distortion_corr_path = subject.get_transform_path(study) + "distortion_corr_32ch"
    try:
        os.makedirs(distortion_corr_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
            

    #    dist_corr settings may differ by study, get them now as needed
    run_settings = study.get_run_type_settings("DISTC")


    #    create NIFTI topup input from concatentated data
    if not os.path.isfile( distortion_corr_path + "/distortion_corr_32ch_APPA.nii.gz" ):

        if not os.path.isfile( distortion_corr_path + "/distortion_corr_32ch_APPA+orig.HEAD" ):
            #    We take the AP PA collection of images and concatenate into a single volume
            #    that will be used to estimate field distortion
            utils.prettyOut("Create distortion_corr_32ch anat volumes...")

            #    AP volume
            if not os.path.isfile( distortion_corr_path + "/distortion_corr_32ch_ap+orig.HEAD" ):
                utils.Cmd("to3d -orient RAI -time:zt "+run_settings['num_slices']+" "+run_settings['time_images']+" "+run_settings['TR']+" seq+z -session " + distortion_corr_path + "/ -prefix distortion_corr_32ch_ap " + path_ap + "/*.dcm").run()

            #    PA volume
            if not os.path.isfile( distortion_corr_path + "/distortion_corr_32ch_pa+orig.HEAD" ):
                utils.Cmd("to3d -orient RAI -time:zt "+run_settings['num_slices']+" "+run_settings['time_images']+" "+run_settings['TR']+" seq+z -session " + distortion_corr_path + "/ -prefix distortion_corr_32ch_pa " + path_pa + "/*.dcm").run()

            #    concatenate to APPA afni
            utils.Cmd("3dTcat -prefix " + distortion_corr_path + "/distortion_corr_32ch_APPA  " + distortion_corr_path + "/distortion_corr_32ch_ap+orig  " + distortion_corr_path + "/distortion_corr_32ch_pa+orig" ).run()

            utils.purge_paths(distortion_corr_path + "/distortion_corr_32ch_ap+orig.*", distortion_corr_path + "/distortion_corr_32ch_pa+orig.*")


        #    convert separtely and carefully to NIFTI format. Never underestimate how inconsistently this is handled by separate AFNI tools!
        utils.Cmd("3dresample -orient RPI -inset " + distortion_corr_path + "/distortion_corr_32ch_APPA+orig -prefix " + distortion_corr_path + "/distortion_corr_32ch_APPA.nii.gz" ).run()



    #    now estimate the field distortion
    if not os.path.isfile( distortion_corr_path + "/dewarped_seEPI_APPA.nii" ):

        utils.prettyOut("estimating field distortions in epi...")

        if not os.path.isfile( study.path + "/params/epi_params.txt" ):
            raise utils.CmdProcError("suscep_estimate: Missing " + study.path + "/params/epi_params.txt")


        utils.Cmd("topup --verbose --imain=" + distortion_corr_path + "/distortion_corr_32ch_APPA --datain=" + study.path + "/params/epi_params.txt --config=" + study.path + "/params/b02b0.cnf --out=" + distortion_corr_path + "/topup_results_APPA --fout=" + distortion_corr_path + "/topup_field_APPA --iout=" + distortion_corr_path + "/dewarped_seEPI_APPA" ).run()

        utils.Cmd("gunzip " + distortion_corr_path + "/dewarped_seEPI_APPA.nii.gz" ).run()



    #    this corrects the header of the dewarped sample only so you can compare it to the input for debugging
    #    otherwise, this is unneccesary in production
    fields     = ['pixdim','srow_x','srow_y','srow_z','qform_code','sform_code','quatern_b','quatern_c','quatern_d','qoffset_x','qoffset_y','qoffset_z']
    utils.correct_niftihdr( distortion_corr_path + "/distortion_corr_32ch_APPA.nii.gz", distortion_corr_path + "/dewarped_seEPI_APPA.nii", fields )


    utils.purge_paths(distortion_corr_path + "/distortion_corr_32ch_*+orig.*")





def correct_field_distortion(subject,study,run_type,param_indx=1):

    """
    apply distortion correction to correct distortions in EPI
    subject - subject object
    study   - study object
    run_type    - MMAT
    param_indx - 1,2 = AP; 3,4 = PA (in current params file)
                 the index here is used to specify the acquisition encoding of the input and match it to the epi_params file.
                 For our current setup, it will be '1' to specify that the entire thing was acquired A -> P
    """


    #    locate distortion correction inputs
    distortion_corr_path = subject.get_transform_path(study) + "distortion_corr_32ch"
    if not os.path.isfile( distortion_corr_path + "/dewarped_seEPI_APPA.nii" ):
        raise utils.CmdProcError("Missing suscep_estimate volume")



    # we need the distortion correction volume to be aligned with the SBREF of this run
    distcorr_to_SBRef(subject,study,run_type)

    #    define run volume names and their paths
    in_vol_names     = (study.label + "." + subject.ursi.short + "." + run_type.upper() + "1.SBRef",
                       study.label + "." + subject.ursi.short + "." + run_type.upper() + ".preproc.uncorrected");
    in_vol_paths        = (subject.get_path(study) + run_type.upper() + "1/",
                        subject.get_path(study) + "REG_VOLUMES/");
    out_vol_names    = (study.label + "." + subject.ursi.short + "." + run_type.upper() + "1.SBRef.dist_corr",
                       study.label + "." + subject.ursi.short + "." + run_type.upper() + ".preproc.dist_corr");


    #    process each volume
    for vol_indx in range(len(in_vol_names)):


        in_vol_name        = in_vol_names[vol_indx]
        vol_path            = in_vol_paths[vol_indx]
        out_vol_name        = out_vol_names[vol_indx]


        #    only do this if final output is not there
        if not os.path.isfile( vol_path + out_vol_name + "+orig.HEAD"):


            #    only do this if have the inputs
            if not os.path.isfile( vol_path + in_vol_name + "+orig.HEAD" ):
                raise utils.CmdProcError("Missing to3d input")

            if not os.path.isfile( distortion_corr_path + "/topup_results_APPA_movpar.txt" ):
                raise utils.CmdProcError("Missing distortion estimate input")

            #    create NIFTI file version of input volume
            if not os.path.isfile( vol_path + in_vol_name + ".nii.gz" ):

                utils.prettyOut(in_vol_name + " : Create NIFTI file version of input volume. This may take awhile. Be patient...")
                cmd = utils.Cmd("3dresample -orient RPI -inset " + vol_path + in_vol_name + "+orig -prefix " + vol_path + in_vol_name + ".nii" )
                cmd.run()


            ##    applytopup --imain=rest.raw.vol0000.nii.gz --inindex=1 --topup=topup_results_APPA --datain=epi_params.txt --method=jac --interp=spline --out=rest.raw.vol0000.dist_corr
            if not os.path.isfile( vol_path + out_vol_name + ".nii.gz" ):

                # use topup files that were aligned to this run type
                run_label        = run_type + "1"
                in_path          = subject.get_path(study) + run_label + "/"
                dc_i_field_file  = in_path + "topup_results_APPA_SBRef"

                utils.prettyOut(in_vol_name + " : Applying distortion correction. This may take awhile. Almost there...")
                cmd = utils.Cmd("applytopup --verbose --imain=" + vol_path + in_vol_name + ".nii --inindex=" + str(param_indx) + " --topup=" + dc_i_field_file + " --datain=" + study.path + "/params/epi_params.txt --method=jac --interp=spline --out=" + vol_path + out_vol_name + ".nii" )
                cmd.run()

#                utils.prettyOut(in_vol_name + " : Applying distortion correction. This may take awhile. Almost there...")
#                cmd = utils.Cmd("applytopup --verbose --imain=" + vol_path + in_vol_name + ".nii --inindex=" + str(param_indx) + " --topup=" + distortion_corr_path + "/topup_results_APPA --datain=" + study.path + "/params/epi_params.txt --method=jac --interp=spline --out=" + vol_path + out_vol_name + ".nii" )
#                cmd.run()


            ## N.B. that spline interpolation can introduce negative values even when given input that has all positive values.
            ## Therefore if ones subsequent analysis methods/software assume positive values one should use something like fslmaths my_splinterpolated_images -abs my_positive_splinterpolated_images before proceeding.
            ##    in AA those are filtered out rather than flipping the sign
            #fslmaths rest.raw.vol0000.dist_corr.nii.gz -thr 0 rest.raw.vol0000.dist_corr.filt -odt short
            if not os.path.isfile( vol_path + out_vol_name + ".filt.nii.gz" ):
                utils.prettyOut(out_vol_name + " : Filter interpolation results")
                cmd = utils.Cmd("fslmaths " + vol_path + out_vol_name + ".nii.gz -thr 0 " + vol_path + out_vol_name + ".filt -odt short" )
                cmd.run()


            if not os.path.isfile( vol_path + out_vol_name + "+orig.HEAD" ):
                utils.prettyOut( out_vol_name  + " : Convert corrected volume to afni file format")
                cmd = utils.Cmd("3dresample -orient RAI -inset " + vol_path + out_vol_name + ".filt.nii.gz -prefix " + vol_path + out_vol_name )
                cmd.run()

            #    movement into NIFTI drops obliquity and sets type to anat, lets correct type here just to avoid warnings when data is viewed later in AFNI
            if os.path.isfile( vol_path + out_vol_name + "+orig.HEAD" ):
                cmd = utils.Cmd("3drefit -epan " + vol_path + out_vol_name + "+orig.HEAD" )
                cmd.run()


            #    clean up
            utils.purge_paths(vol_path + in_vol_name + "*nii*",
                              vol_path + out_vol_name + "*nii*",
                              subject.get_path(study) + "REG_VOLUMES/" + study.label + "." + subject.ursi.short + "." + run_type.upper() + ".preproc.uncorrected.*")







def qa_correct_field_distortion(subject,study,run_type):

    """
    QA display of field distortion

    """

    #utils.prettyOut("epi.py - qa_correct_field_distortion() no longer available as .preproc.uncorrected+orig is no longer available")
    #return

    #    this bit checks the APPA distortion correction out of topup
    ##    locate distortion correction inputs
    #distortion_corr_path = subject.get_transform_path(study) + "distortion_corr_32ch"
    #if not os.path.isfile( distortion_corr_path + "/dewarped_seEPI_APPA.nii" ):
        #raise utils.CmdProcError("Missing suscep_estimate volume")

    #os.chdir(distortion_corr_path)

    ##    Don't let afni session detach or it will close before user can see images
    #os.environ["AFNI_DETACH"] = "NO"

    #cmd = utils.Cmd('afni -com "SWITCH_UNDERLAY distortion_corr_32ch_APPA.nii.gz" -com "SWITCH_OVERLAY dewarped_seEPI_APPA.nii.gz" -com "SEE_OVERLAY +" -com "SET_THRESHOLD A.0" -com "SET_PBAR_ALL A.+99 1.0" -com "SET_SPM_XYZ 6 12 14" ')
    #utils.cmd_echo(cmd)
    #return_code = subprocess.call(cmd, shell=True)
    #if return_code:
        #raise utils.CmdProcError(return_code)

    ##    return this setting to lab default
    #os.environ["AFNI_DETACH"] = "YES"



    os.chdir(subject.get_path(study) + "REG_VOLUMES/")

    #    full volume load time is slow. Don't know how to display in AFNI without loading entire volume so i'll create temp volumes for inspection
    ul_file_path = os.path.join(subject.get_path(study), run_type.upper() + "1/", study.label + "." + subject.ursi.short + "." + run_type.upper() + "1.SBRef+orig")
    utils.Cmd('3dcalc -a ' + ul_file_path + ' -prefix tmp.verify.underlay -expr a').run()
    utils.Cmd('3dcalc -a ' + study.label + '.' + subject.ursi.short + '.' + run_type + '.preproc.dist_corr+orig"[0]"   -prefix tmp.verify.overlay  -expr a').run()

    #    Don't let afni session detach or it will close before user can see images
    os.environ["AFNI_DETACH"] = "NO"

    utils.Cmd('afni -com "SET_ANATOMY tmp.verify.underlay " -com "SET_FUNCTION tmp.verify.overlay" -com "SEE_OVERLAY +" -com "SET_THRESHOLD A.0" -com "SET_PBAR_ALL A.+99 1.0" -com "SET_SPM_XYZ 6 12 14" ').run()

    #    return this setting to lab default
    os.environ["AFNI_DETACH"] = "YES"

    utils.purge_paths("tmp.verify.*lay*")







def create_run_motion_regressors(subject,study,run_type,base=1):

    """
    creates file appropriate for level1

    """

    import glob

    out_path = study.movement_path + "/" + subject.visit
    try:
        os.makedirs(out_path,0o770)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    #    path for this run type
    run_path = subject.get_path(study) + run_type
    #    number of runs for this run type
    run_count = study.get_run_type_settings(run_type)['run_count']


    #     original movement regressor created with raw epi input. This is used for motion quantification
    #    create from available run data

    if not os.path.isfile(out_path + "/" + subject.ursi.short + "." + run_type + ".Movement.Regressor.1D"):
        log_files = sorted(glob.glob(run_path + "?/*."+run_type+"*.3dreglog.raw"))

        #    errors for nshaff on tethys which might be caused by written file not yet being ready for read
        #    to debug, count and stop. nshaff should then rerun to see if that works
        print(("Run Count: {}".format(run_count)))
        if (len(log_files) != int(run_count)):
            raise utils.CmdProcError("Not enough 3dreglog.raw input files.")

        for log in log_files:
            cmd = utils.Cmd("tail -n +" + str(base+1) + " " + log +" >> " + out_path + "/" + subject.ursi.short + "." + run_type + ".Movement.Regressor.1D")
            cmd.run()


    #    additional movement regressor created with corrected epi input. This is used as motion regressor
    if not os.path.isfile(out_path + "/" + subject.ursi.short + "." + run_type + ".Movement.Regressor.mc.1D"):
        log_files = sorted(glob.glob(run_path + "?/*."+run_type+"*.3dreglog.mc"))
        #    errors for nshaff on tethys which might be caused by written file not yet being ready for read
        #    to debug, count and stop. nshaff should then rerun to see if that works
        if (len(log_files) != int(run_count)):
            raise utils.CmdProcError("Not enough 3dreglog.mc input files.")

        for log in log_files:
            cmd = utils.Cmd("tail -n +" + str(base+1) + " " + log +" >> " + out_path + "/" + subject.ursi.short + "." + run_type + ".Movement.Regressor.mc.1D")
            cmd.run()


    #    create derivative version, aka frame-wise displacement
    #    regression requires derivative of motion regressor
    #    this has been failing because above concatenation not happening correctly.
    if not os.path.isfile(out_path + "/" + subject.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D"):

        run_count = study.get_run_type_settings(run_type)['run_count']

        #    sadly, 1d_tool doesn't appear to be importable. i call it through system
        cmd = utils.Cmd("1d_tool.py -infile " + out_path + "/" + subject.ursi.short + "." + run_type + ".Movement.Regressor.mc.1D -set_nruns " + str(run_count) + " -derivative -write " + out_path + "/" + subject.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D")
        cmd.run()






def clean_run_dir(subj,study,run_type):

    """
    clean up stuff we don't need after final preprocessed volume

    """
    # intermediate step output is removed here at the end
    # the '' empty pattern removes the original volume converted from dicom
    for match_part in ['.dspk', '.tshft', '.2dreg', '.3dreg', '']:
        utils.prettyOut(subj.get_path(study) + run_type + "?/" + study.label + "." + subj.ursi.short + "." + run_type + "?" + match_part + "+orig*")
        utils.purge_paths(subj.get_path(study) + run_type + "?/" + study.label + "." + subj.ursi.short + "." + run_type + "?" + match_part + "+orig*")

    # others
    utils.purge_paths(subj.get_path(study) + run_type + "1/" + study.label + "." + subj.ursi.short + "." + run_type + "1.SBRef.nii.gz")







def generateResiduals(subj,study,fwhm,run_label,base_suffix="dist_corr"):

    """
    generate residuals for resting state analyzes
    Deprecate?
    rsfc_create_residuals() can generate residuals with or without blurring prior
    The only difference in this function then is that bandpass filters and spatially normalizes the residuals

    """
    #    build subject_path
    subject_path = study.path + "/" + subj.ursi.full

    #    workaround fact that regressor files contain visit
    #    and any other differences between studies that have or have not a visit....
    #subject_visit_path = ""
    if subj.visit:
        #subject_visit_path = "." + subj.visit
        subject_path = subject_path + "/" + subj.visit


    #    run_label indicates where processing data is
    #    run_type indicates final processing outpute
    run_type = run_label.rstrip("0123456789")
    final_blur_fwhm = str(fwhm)


    #    check for final output only before running
    if not os.path.isfile(subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3+tlrc.HEAD"):

        #    deposit junk in /tmp
        os.chdir("/tmp")

        #    define inputs to check

        # allow null base suffix
        if len(base_suffix):
            base_suffix = "."+base_suffix
        rest_input_file = subject_path + "/REG_VOLUMES/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc"+base_suffix+"+orig"
        rest_tx_epi2T1_file = subj.get_transform_path(study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2T1.aff12.1D"
        rest_tx_T12epi_file = subj.get_transform_path(study) + study.label + "." + subj.ursi.short + "." + run_type + "_T12epi.aff12.1D"
        rest_movement_file = study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.1D"
        spm_wm_file = subject_path + "/anat/" + study.label +"."+ subj.ursi.short+".anat.spm.wm.mask.60+orig"
        spm_csf_file = subject_path + "/anat/" + study.label +"."+ subj.ursi.short+".anat.spm.csf.mask.60+orig"



        #    check inputs
        if not os.path.isfile(rest_input_file+".HEAD"):

            raise utils.CmdProcError("\nMissing Input " + rest_input_file + "\n")


        if not os.path.isfile(rest_tx_T12epi_file):

            if not os.path.isfile(rest_tx_epi2T1_file):

                raise utils.CmdProcError("\nMissing Input " + rest_tx_epi2T1_file + "\n")

            else:

                cmd = utils.Cmd("cat_matvec -ONELINE " + rest_tx_epi2T1_file + " -I > " + rest_tx_T12epi_file )
                cmd.run()



        if not os.path.isfile(rest_movement_file):

            raise utils.CmdProcError("\nMissing Input " + rest_movement_file + "\n")


        if not os.path.isfile(spm_wm_file+".HEAD"):

            raise utils.CmdProcError("\nMissing Input " + spm_wm_file + "\n")


        if not os.path.isfile(spm_csf_file+".HEAD"):

            raise utils.CmdProcError("\nMissing Input " + spm_csf_file + "\n")






        #    now start creating the tissue regressor files
        utils.prettyOut("Resting State Residuals : Creating the tissue regressor files")
        if not os.path.isfile( subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.wm.1D" ):


            if not os.path.isfile( subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.wm." + run_type + "_Reg+orig.HEAD" ):

                cmd = utils.Cmd("3dAllineate "
                        "-master " + rest_input_file + " "
                        "-1Dmatrix_apply " + rest_tx_T12epi_file + " "
                        "-input " + spm_wm_file +" "
                        "-final NN "
                        "-prefix "+subject_path + "/anat/" + study.label +"."+ subj.ursi.short+".anat.spm.wm." + run_type + "_Reg" )

                cmd.run()

            cmd = utils.Cmd("3dROIstats -quiet "
                    "-mask " +subject_path + "/anat/" + study.label + "." + subj.ursi.short+ ".anat.spm.wm." + run_type + "_Reg+orig "
                    + rest_input_file + " > " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.wm.1D" )
            cmd.run()



        if not os.path.isfile( subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.csf.1D" ):

            if not os.path.isfile( subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.spm.csf." + run_type + "_Reg+orig.HEAD" ):
                cmd = utils.Cmd("3dAllineate "
                        "-master " + rest_input_file + " "
                        "-1Dmatrix_apply " + rest_tx_T12epi_file + " "
                        "-input " + spm_csf_file +" "
                        "-final NN "
                        "-prefix "+subject_path + "/anat/" + study.label +"."+ subj.ursi.short+".anat.spm.csf." + run_type + "_Reg" )
                cmd.run()


            cmd = utils.Cmd("3dROIstats -quiet "
                    "-mask " +subject_path + "/anat/" + study.label + "." + subj.ursi.short+ ".anat.spm.csf." + run_type + "_Reg+orig "
                    + rest_input_file + " > " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.csf.1D" )
            cmd.run()

            #cmd = utils.Cmd("rm -f " +subject_path + "/anat/" + study.label + "." + subj.ursi.short+ ".anat.spm.csf." + run_type + "_Reg+orig* " )
            #return_code = subprocess.call(cmd, shell=True)
            #if return_code:
                #raise utils.CmdProcError(return_code)



        #    generate motion derivative regessors if they don't exist
        if not os.path.isfile(study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D"):

            cmd = utils.Cmd("1d_tool.py"
                    " -infile " + rest_movement_file +
                    " -set_nruns 1"
                    " -derivative"
                    " -write " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D" )

            cmd.run()



        #    blur input to glm
        if not os.path.isfile( subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm + "+orig.HEAD" ):

            cmd = utils.Cmd("3dmerge"
                " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm +
                " -1blur_fwhm "  + final_blur_fwhm +
                " -doall "
                + rest_input_file )

            cmd.run()



        #    run glm to generate residuals
        utils.prettyOut("Resting State Residuals : Run glm to generate residuals")
        if not os.path.isfile( subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + "+orig.HEAD" ):

            cmd = utils.Cmd("3dDeconvolve"
                    " -overwrite"
                    " -input " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm + "+orig"
                    " -polort 2"
                    " -allzero_OK"
                    " -GOFORIT 10"
                    " -num_stimts 14"
                    " -stim_file 1 " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.csf.1D"
                    " -stim_label 1 'CSF'"
                    " -stim_minlag 1 0"
                    " -stim_maxlag 1 0"
                    " -stim_file 2 " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + ".Tissue.Regressor.wm.1D"
                    " -stim_label 2 'WM'"
                    " -stim_minlag 2 0"
                    " -stim_maxlag 2 0"
                    " -stim_file 3 " + rest_movement_file +"'[1]'"
                    " -stim_label 3 'ROLL'"
                    " -stim_minlag 3 0"
                    " -stim_maxlag 3 0"
                    " -stim_file 4 " + rest_movement_file +"'[2]'"
                    " -stim_label 4 'PITCH'"
                    " -stim_minlag 4 0"
                    " -stim_maxlag 4 0"
                    " -stim_file 5 " + rest_movement_file +"'[3]'"
                    " -stim_label 5 'YAW'"
                    " -stim_minlag 5 0"
                    " -stim_maxlag 5 0"
                    " -stim_file 6 " + rest_movement_file +"'[4]'"
                    " -stim_label 6 'dS'"
                    " -stim_minlag 6 0"
                    " -stim_maxlag 6 0"
                    " -stim_file 7 " + rest_movement_file +"'[5]'"
                    " -stim_label 7 'dL'"
                    " -stim_minlag 7 0"
                    " -stim_maxlag 7 0"
                    " -stim_file 8 " + rest_movement_file +"'[6]'"
                    " -stim_label 8 'dP'"
                    " -stim_minlag 8 0"
                    " -stim_maxlag 8 0"
                    " -stim_file 9 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[1]'"
                    " -stim_label 9 'ROLL_deriv'"
                    " -stim_minlag 9 0"
                    " -stim_maxlag 9 0"
                    " -stim_file 10 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[2]'"
                    " -stim_label 10 'PITCH'_deriv"
                    " -stim_minlag 10 0"
                    " -stim_maxlag 10 0"
                    " -stim_file 11 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[3]'"
                    " -stim_label 11 'YAW_deriv'"
                    " -stim_minlag 11 0"
                    " -stim_maxlag 11 0"
                    " -stim_file 12 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[4]'"
                    " -stim_label 12 'dS_deriv'"
                    " -stim_minlag 12 0"
                    " -stim_maxlag 12 0"
                    " -stim_file 13 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[5]'"
                    " -stim_label 13 'dL_deriv'"
                    " -stim_minlag 13 0"
                    " -stim_maxlag 13 0"
                    " -stim_file 14 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[6]'"
                    " -stim_label 14 'dP_deriv'"
                    " -stim_minlag 14 0"
                    " -stim_maxlag 14 0"
                    " -errts " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm +
                    " -jobs 10 "
                    " -xsave"
                    " -float"
                    " -nobucket" )

            cmd.run()



        #    bandpass filtering
        utils.prettyOut("Resting State Residuals : bandpass filtering")
        if not os.path.isfile( subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1+orig.HEAD" ):

            cmd = utils.Cmd(" 3dBandpass"
                    " -band 0.01 0.1"
                    " -input " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + "+orig"
                    " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1" )

            cmd.run()



        # convert residual file to tlrc at 3x3x3 resolution
        utils.prettyOut("Resting State Residuals : convert residual file to tlrc at 3x3x3 resolution")
        if not os.path.isfile( subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3+tlrc.HEAD" ):

            #cmd = utils.Cmd(" 3dAllineate"
                    #" -overwrite"
                    #" -master " + subject_path + "/anat/" + study.label + "." + subj.ursi.short + ".anat.sksp+tlrc"
                    #" -mast_dxyz 3.0"
                    #" -1Dmatrix_apply " + subj.get_transform_path(study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2tlrc.aff12.1D"
                    #" -input " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1+orig"
                    #" -floatize"
                    #" -final trilinear"
                    #" -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3" )


            nl_warp_src        = subj.get_transform_path(study) + study.label+"."+subj.ursi.short+".anat.sksp.qwarp_WARP1+tlrc"
            t1_to_tlrc_mat    = subj.get_transform_path(study) + study.label+"."+subj.ursi.short+".T1_2tlrc.aff3x4.1D"
            epi_2_T1_aff_mat    = subj.get_transform_path(study) + study.label+"."+subj.ursi.short+"."+run_type+"_epi2T1.aff12.1D"

            #    the reference volume output will be aligned to
            #    anat.sksp.qwarp+tlrc
            #    which is non-linear warp to TT_N27
            #    so reference volume output should look good over that
            cmd = utils.Cmd("3dNwarpApply " +
                    " -source " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1+orig " +
                    " -nwarp '" + nl_warp_src + " " + t1_to_tlrc_mat + " " + epi_2_T1_aff_mat + "'" +
                    " -master " + subj.get_path(study) + "anat/" + study.label+"."+subj.ursi.short+".anat.sksp.qwarp+tlrc" +
                    " -interp wsinc5" +
                    " -dxyz 3.0" +
                    " -prefix " + subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".resid.fwhm" + final_blur_fwhm + ".bp-0_01-0_1.3x3x3" )

            cmd.run()



    #    clean up even if nothing else was run
    utils.purge_paths(subject_path + "/" + run_label + "/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc.fwhm" + final_blur_fwhm + "+orig*",
                    subject_path + "/anat/" + study.label + "." + subj.ursi.short+ ".anat.spm.*." + run_type + "_Reg+orig* ")







def create_tissue_regressor(subj, study, run_label, tissue_classes=['wm', 'csf'], base_suffix="dist_corr"):

    """
    Generate tissue regressor files for any run type.
    Regressor file outputs into run_label directory

    """

    #    run_label indicates where processing data is
    #    run_type indicates final processing output
    run_type = run_label.rstrip("0123456789")


    #    define inputs to check

    # allow null base suffix
    if len(base_suffix):
        base_suffix = "."+base_suffix
        
    
    rest_input_file = subj.get_path(study) + "/REG_VOLUMES/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc"+base_suffix+"+orig.HEAD"
    rest_tx_T12epi_file = subj.get_transform_path(study) + study.label + "." + subj.ursi.short + "." + run_type + "_T12epi.aff12.1D"
    rest_tx_epi2T1_file = subj.get_transform_path(study) + study.label + "." + subj.ursi.short + "." + run_type + "_epi2T1.aff12.1D"
    
    

    for tissue_class in tissue_classes:
        
        spm_tissue_file = os.path.join( subj.get_path(study, 'anat'), subj.id_fn(study,fn="anat.spm."+tissue_class+".mask.60+orig.HEAD") )
        if not os.path.isfile(spm_tissue_file):
            raise utils.CmdProcError("\nMissing Input " + spm_tissue_file + "\n")            

        #    check preproc input
        if not os.path.isfile(rest_input_file):
            raise utils.CmdProcError("\nMissing Input " + rest_input_file + "\n")
            
        #   check that transform exists
        if not os.path.isfile(rest_tx_T12epi_file):
            if not os.path.isfile(rest_tx_epi2T1_file):
                raise utils.CmdProcError("\nMissing Input " + rest_tx_epi2T1_file + "\n")
            else:
                utils.Cmd("cat_matvec -ONELINE " + rest_tx_epi2T1_file + " -I > " + rest_tx_T12epi_file ).run()


        #    now start creating the tissue regressor files
        o_file_1d = os.path.join( subj.get_path(study), run_label, subj.id_fn(study,fn="Tissue.Regressor."+tissue_class+".1D") )
        o_file_reg = os.path.join( subj.get_path(study, 'anat'), subj.id_fn(study,fn="anat.spm." + tissue_class + "." + run_type + "_Reg+orig.HEAD") )
                
        if not os.path.isfile( o_file_1d ):            
            if not os.path.isfile( o_file_reg ):
                utils.prettyOut("RSFC : Creating the tissue regressor files")
                utils.Cmd("3dAllineate "
                        "-master " + rest_input_file + " "
                        "-1Dmatrix_apply " + rest_tx_T12epi_file + " "
                        "-input " + spm_tissue_file +" "
                        "-final NN "
                        "-prefix "+o_file_reg ).run()


            utils.Cmd("3dROIstats -quiet -mask " + o_file_reg + " " + rest_input_file + " > " + o_file_1d ).run()

            utils.purge_paths(o_file_reg)
            





def rsfc_create_residuals(subj, study, run_label, fwhm=0, base_suffix="dist_corr"):

    """
    Generate residuals. Fwhm is optional

    """

    #    run_label indicates where processing data is
    #    run_type indicates final processing outpute
    run_type = run_label.rstrip("0123456789")
    final_blur_fwhm = str(fwhm)


    final_o_file = os.path.join(subj.get_path(study), run_label, subj.id_fn(study,fn=run_type + ".resid+orig.HEAD") )
    if fwhm != 0:
        final_o_file = final_o_file.replace(".resid", ".resid.fwhm" + final_blur_fwhm)
    
    
    #    check for final output only before running
    if not os.path.isfile( final_o_file + ".HEAD"):

        # allow null base suffix
        if len(base_suffix):
            base_suffix = "."+base_suffix
        
        #    define inputs to check        
        i_file = subj.get_path(study) + "/REG_VOLUMES/" + study.label + "." + subj.ursi.short + "." + run_type + ".preproc"+base_suffix+"+orig.HEAD"
        rest_movement_file = study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.1D"
        #    check preproc input
        if not os.path.isfile(i_file):
            raise utils.CmdProcError("\nMissing Input " + i_file + "\n")
        if not os.path.isfile(rest_movement_file):
            raise utils.CmdProcError("\nMissing Input " + rest_movement_file + "\n")



        # make sure we have tissue regressors
        create_tissue_regressor(subj, study, run_label)
        # make sure we have motion derivatives
        create_motion_deriv_file(subj, study, run_label)


        # blur input?
        if fwhm != 0:
            i_file_fwhm = i_file.replace(".preproc", ".preproc.fwhm" + final_blur_fwhm)
            
            if not os.path.isfile( i_file_fwhm ):
                utils.Cmd("3dmerge"
                    " -prefix " + i_file_fwhm + 
                    " -1blur_fwhm "  + final_blur_fwhm +
                    " -doall "
                    + i_file ).run()
    
            i_file = i_file_fwhm



        #    run glm to generate residuals
        if not os.path.isfile( final_o_file ):

            utils.prettyOut("Resting State Residuals : Run glm to generate residuals")
            
            #    deposit junk in /tmp
            os.chdir("/tmp")
            
            utils.Cmd("3dDeconvolve"
                    " -overwrite"
                    " -input " + i_file +
                    " -polort 2"
                    " -allzero_OK"
                    " -GOFORIT 10"
                    " -num_stimts 14"
                    " -stim_file 1 " + os.path.join( subj.get_path(study), run_label, subj.id_fn(study,fn="Tissue.Regressor.csf.1D") ) +
                    " -stim_label 1 'CSF'"
                    " -stim_minlag 1 0"
                    " -stim_maxlag 1 0"
                    " -stim_file 2 " + os.path.join( subj.get_path(study), run_label, subj.id_fn(study,fn="Tissue.Regressor.wm.1D") ) +
                    " -stim_label 2 'WM'"
                    " -stim_minlag 2 0"
                    " -stim_maxlag 2 0"
                    " -stim_file 3 " + rest_movement_file +"'[1]'"
                    " -stim_label 3 'ROLL'"
                    " -stim_minlag 3 0"
                    " -stim_maxlag 3 0"
                    " -stim_file 4 " + rest_movement_file +"'[2]'"
                    " -stim_label 4 'PITCH'"
                    " -stim_minlag 4 0"
                    " -stim_maxlag 4 0"
                    " -stim_file 5 " + rest_movement_file +"'[3]'"
                    " -stim_label 5 'YAW'"
                    " -stim_minlag 5 0"
                    " -stim_maxlag 5 0"
                    " -stim_file 6 " + rest_movement_file +"'[4]'"
                    " -stim_label 6 'dS'"
                    " -stim_minlag 6 0"
                    " -stim_maxlag 6 0"
                    " -stim_file 7 " + rest_movement_file +"'[5]'"
                    " -stim_label 7 'dL'"
                    " -stim_minlag 7 0"
                    " -stim_maxlag 7 0"
                    " -stim_file 8 " + rest_movement_file +"'[6]'"
                    " -stim_label 8 'dP'"
                    " -stim_minlag 8 0"
                    " -stim_maxlag 8 0"
                    " -stim_file 9 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[1]'"
                    " -stim_label 9 'ROLL_deriv'"
                    " -stim_minlag 9 0"
                    " -stim_maxlag 9 0"
                    " -stim_file 10 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[2]'"
                    " -stim_label 10 'PITCH'_deriv"
                    " -stim_minlag 10 0"
                    " -stim_maxlag 10 0"
                    " -stim_file 11 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[3]'"
                    " -stim_label 11 'YAW_deriv'"
                    " -stim_minlag 11 0"
                    " -stim_maxlag 11 0"
                    " -stim_file 12 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[4]'"
                    " -stim_label 12 'dS_deriv'"
                    " -stim_minlag 12 0"
                    " -stim_maxlag 12 0"
                    " -stim_file 13 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[5]'"
                    " -stim_label 13 'dL_deriv'"
                    " -stim_minlag 13 0"
                    " -stim_maxlag 13 0"
                    " -stim_file 14 " + study.movement_path + "/" + subj.visit + "/" + subj.ursi.short + "." + run_type + ".Movement.Regressor.mc.Deriv.1D'[6]'"
                    " -stim_label 14 'dP_deriv'"
                    " -stim_minlag 14 0"
                    " -stim_maxlag 14 0"
                    " -errts " + final_o_file +
                    " -jobs 6 "
                    " -xsave"
                    " -float"
                    " -nobucket" ).run()


        if fwhm != 0:
            #    clean up even if nothing else was run
            utils.purge_paths(i_file)






def rsfc_reho_falff(subj, study, fwhm):

    """
    Regional homogeniety (ReHo)
    Fractional amplitude of low frequency fluctuations (falff) : Similar to alff except that the variance of the bandpassed signal is 
    divided by the total variance (variance of non-bandpassed signal).
    fALFF is a ratio:
        (the sum of amplitudes of those LFFs) / (sum of amplitudes of all frequencies in time series). 

    See Zou 2008
    """
    import subprocess
    
    #    build subject_path
    subject_path = subj.get_path(study)
    final_blur_fwhm = str(fwhm)


    i_mask = os.path.join(subject_path, "REG_VOLUMES", study.label + "." + subj.ursi.short + ".REST.BASEIMG_sksp+orig.HEAD")
    o_mask = os.path.join(subject_path, "REG_VOLUMES", study.label + "." + subj.ursi.short + ".REST.BASEIMG_mask_dil_1+orig.HEAD")
    if not os.path.isfile(o_mask):
        utils.Cmd("3dmask_tool -input "+i_mask+" -prefix "+o_mask+" -dilate_input 1").run()
        
    
    #    check for final output only before running
    i_file = os.path.join(subject_path, 'REST1', study.label + "." + subj.ursi.short + ".REST.resid+orig.HEAD" )
    o_file = os.path.join(subject_path, 'REST1', study.label + "." + subj.ursi.short + ".REST.resid.rsfc" )
    if not os.path.isfile(o_file + '_LFF+orig.HEAD'):
        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : bandpass and calculate metrics" )
        utils.Cmd("3dRSFC -prefix " + o_file + " -mask " + o_mask + " -nodetrend -band 0.01 0.1 -input " + i_file ).run()
        
    
    i_file = os.path.join(subject_path, 'REST1', study.label + "." + subj.ursi.short + ".REST.resid.rsfc_LFF+orig.HEAD" )
    o_file = os.path.join(subject_path, 'REST1', study.label + "." + subj.ursi.short + ".REST.resid.rsfc_REHO+orig.HEAD"  )
    if not os.path.isfile(o_file):
        utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : calculate REHO on bandpass filtered data" )
        utils.Cmd("3dReHo -prefix " + o_file + " -mask " + o_mask + "  -inset " + i_file ).run()
        
        
    
    
    # blur and normalize metrics of interest
    for metric in ['REHO', 'fALFF']:
        
        final_o_metric = os.path.join(subject_path, 'REST1', study.label + "." + subj.ursi.short + ".REST.resid.rsfc_" + metric + "_Z.fwhm" + final_blur_fwhm + ".1x1x1+tlrc.HEAD" )
        if not os.path.isfile(final_o_metric):
            
            # now blur
            if fwhm != 0:
                i_file = os.path.join(subject_path, 'REST1', study.label + "." + subj.ursi.short + ".REST.resid.rsfc_" + metric + "+orig.HEAD" )
                o_file = i_file.replace(".rsfc_"+metric, ".rsfc_"+metric+".fwhm"+final_blur_fwhm)
    
                if not os.path.isfile( o_file ):
                    utils.Cmd("3dmerge"
                        " -prefix " + o_file + 
                        " -1blur_fwhm "  + final_blur_fwhm +
                        " -doall "
                        + i_file ).run()
        
                i_file = o_file
            
            
            # convert to z. See Zou 2008
            o_file = os.path.join(subject_path, 'REST1', study.label + "." + subj.ursi.short + ".REST.resid.rsfc_" + metric + "_Z" + ".fwhm" + final_blur_fwhm + "+orig.HEAD" )
            if not os.path.isfile(o_file):
                utils.prettyOut(subj.ursi.full + " : " + subj.visit + " : " + metric + " : Convert to Z-scores" )
                #utils.Cmd("3dmaskave -quiet -sigma -mask "+ mask + " " + i_file).echo()
                cmd_str = "3dmaskave -quiet -sigma -mask "+ o_mask + " " + i_file
                subproc_res = subprocess.run(cmd_str, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                mean = subproc_res.stdout.split('\n')[1].split()[0]
                sd = subproc_res.stdout.split('\n')[1].split()[1]
                
                utils.Cmd("3dcalc -a "+i_file+" -b "+o_mask+" -expr 'b*((a-"+mean+")/"+sd+")' -prefix " + o_file).run()
                
                utils.purge_paths(i_file.replace(".HEAD","*"))
            
    
            # tlrc at 1x1x1 resolution
            i_file = o_file    
            nl_warp_src        = subj.get_transform_path(study) + study.label+"."+subj.ursi.short+".anat.sksp.qwarp_WARP1+tlrc"
            t1_to_tlrc_mat     = subj.get_transform_path(study) + study.label+"."+subj.ursi.short+".T1_2tlrc.aff3x4.1D"
            epi_2_T1_aff_mat   = subj.get_transform_path(study) + study.label+"."+subj.ursi.short+".REST_epi2T1.aff12.1D"

            #    the reference volume output will be aligned to
            #    anat.sksp.qwarp+tlrc
            #    which is non-linear warp to TT_N27
            #    so reference volume output should look good over that
            utils.Cmd("3dNwarpApply " +
                    " -source " + i_file +
                    " -nwarp '" + nl_warp_src + " " + t1_to_tlrc_mat + " " + epi_2_T1_aff_mat + "'" +
                    " -master " + subj.get_path(study) + "anat/" + study.label+"."+subj.ursi.short+".anat.sksp.qwarp+tlrc" +
                    " -interp wsinc5" +
                    " -prefix " + final_o_metric ).run()

            utils.purge_paths(i_file.replace(".HEAD","*"))
            



