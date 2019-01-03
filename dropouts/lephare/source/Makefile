#----------------------------------------------------------
#  ZPHOT
#    Makefile
#----------------------------------------------------------
#
#  UNIT
#
#       File:	Makefile 
#
#       Author:	Benoit Vandame    :  Many Thanks guy !!
#
#       Date:	06/03/2000
#     
#----------------------------------------------------------
#
#  DESCRIPTION
#
#       Makefile for LE PHARE  Programs
#
#----------------------------------------------------------
#
#
################################################
#######       Compilation      #################
################################################
#
#
#######    For Other UNIX boxes      ###########
#FC         = f77 
#FFLAGS     = -O5 
#LIBS       = 
#
#######         For Linux              #########
FC         = gfortran
#FC         = g77
FFLAGS     = -O5 -Wall 
LIBS       = 
#LIBS       = -lg2c -lm
#
#######  For debugging (Linux/Solaris) ##########
#FFLAGS     =  -pg -g  
#
########       Create the LIBRARY      ##########
LIB_LEPHARE  = lephare_lib.a
AR         = ar
RANLIB     = ranlib
#
########         Source archive        ##########
ARCHIVE = lephare-1.0.0.tar
#
########        Compiling source       ##########

FILES  = \
	sedtolib\
	filter\
	filter_info\
	filter_extinc\
	mag_star\
        mag_gal\
        zphota\
	mag_zform
#        zphota_v0\
#	mag_eros\
#        zphota_sc\
#        zphota_fir\
#	dz_color\
#	simul\
#        read_specz0
       
OBJECTS = \
	prep_output.o\
	check_libphys.o\
	chi_libphys.o\
	chi_para.o\
	nzpriorVVDStest07.o\
        nzpriorAGN07.o\
        ext_emlines.o\
        fluxem.o\
        addEmission.o\
        addEmission2.o\
	addemlines.o\
        fit_PDFz.o\
	read_para.o\
	read_para2.o\
	get_conf.o\
	getc_option.o\
	geti_option.o\
	getf_option.o\
	get_value.o\
	get_path.o\
	opt_line.o\
	val_string.o\
	check_float.o\
	help.o\
	get_help.o\
	err_option.o\
	sampling.o\
	lbref.o\
	lambda.o\
	scale_ext.o\
	scale_opa.o\
	indexx.o\
	indexz.o\
	cal_mag.o\
	funz.o\
	time.o\
	dmet.o\
	zage.o\
	zeropoint.o\
	context.o\
	probaz.o\
	trapzd.o\
	probazbay.o\
	dchi2.o\
	zpeak_scale.o\
	int_parab.o\
	sort_color.o\
	search_rec.o\
	locate.o\
	gammln.o\
        auto_adapt.o\
	sortNR.o\
        fcnfunc.o\
	nzprior.o\
	nzprior2.o\
	nzprior3.o\
	nzprior4.o\
	nzprior5.o\
	getOpt1Tab.o\
	getOpt2Tab.o\
	abs_mag.o\
	abs_mag_onefilt.o\
	k_interp.o\
	z_vmax.o\
        xweight.o\
	read_lib.o\
	read_libfir.o\
	read_sclib.o\
        chi_sc.o\
        read_spec.o\
        read_star.o\
	write_out.o\
	write_sim.o\
        minuit/mnamin.o\
        minuit/mnbins.o\
        minuit/mncalf.o\
        minuit/mncler.o\
        minuit/mncntr.o\
        minuit/mncomd.o\
        minuit/mncont.o\
        minuit/mncrck.o\
        minuit/mncros.o\
        minuit/mncuve.o\
        minuit/mnderi.o\
        minuit/mndxdi.o\
        minuit/mneig.o\
        minuit/mnemat.o\
        minuit/mnerrs.o\
        minuit/mneval.o\
        minuit/mnexcm.o\
        minuit/mnexin.o\
        minuit/mnfixp.o\
        minuit/mnfree.o\
        minuit/mngrad.o\
        minuit/mnhelp.o\
        minuit/mnhes1.o\
        minuit/mnhess.o\
        minuit/mnimpr.o\
        minuit/mninex.o\
        minuit/mninit.o\
        minuit/mninpu.o\
        minuit/mnintr.o\
        minuit/mnlims.o\
        minuit/mnline.o\
        minuit/mnmatu.o\
        minuit/mnmigr.o\
        minuit/mnmnos.o\
        minuit/mnmnot.o\
        minuit/mnparm.o\
        minuit/mnpars.o\
        minuit/mnpfit.o\
        minuit/mnpint.o\
        minuit/mnplot.o\
        minuit/mnpout.o\
        minuit/mnprin.o\
        minuit/mnpsdf.o\
        minuit/mnrazz.o\
        minuit/mnread.o\
        minuit/mnrn15.o\
        minuit/mnrset.o\
        minuit/mnsave.o\
        minuit/mnscan.o\
        minuit/mnseek.o\
        minuit/mnset.o\
        minuit/mnseti.o\
        minuit/mnsimp.o\
        minuit/mnstat.o\
        minuit/mnstin.o\
        minuit/mntiny.o\
        minuit/mnunpt.o\
        minuit/mnvers.o\
        minuit/mnvert.o\
        minuit/mnwarn.o\
        minuit/mnwerr.o\
        minuit/stand.o
all: $(FILES) work

%.o: %.f
	$(FC) $(FFLAGS) -c  $*.f  -o $*.o 

%: %.f $(LIB_LEPHARE)
	$(FC) $(FFLAGS) $*.f -o $*  $(LIB_LEPHARE) $(LIBS) 

$(FILES): $(LIB_LEPHARE)

$(LIB_LEPHARE): $(OBJECTS)
	rm -f $(LIB_LEPHARE)
	$(AR) rv $(LIB_LEPHARE) $(OBJECTS)
	$(RANLIB) $(LIB_LEPHARE)
#
#
#################################################
########            Cleaning             ########
#################################################
#
#
clean:
	rm -f \#*
	rm -f *.o
	rm -f *.a
	rm -f *~

	rm -f minuit/\#*
	rm -f minuit/*.o
	rm -f minuit/*.a
	rm -f minuit/*~

#
#
#################################################
########            Archiving           #########
#################################################
#
#
archi:
	rm -f $(ARCHIVE)
	tar cvf $(ARCHIVE) *.f
	gzip $(ARCHIVE)
	echo "Archive: $(ARCHIVE).gz"
#
#
################################################
#######     Create working directory        ####
################################################
work:
	mkdir -p $(LEPHAREWORK)/filt
	mkdir -p $(LEPHAREWORK)/lib_bin
	mkdir -p $(LEPHAREWORK)/lib_mag
#
