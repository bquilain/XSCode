#-- start of make_header -----------------

#====================================
#  Library MyIngridLibs
#
#   Generated Wed Feb  6 18:49:03 2013  by bquilain
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_MyIngridLibs_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_MyIngridLibs_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_MyIngridLibs

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_MyIngridLibs = /tmp/CMT_$(INGRID_tag)_MyIngridLibs.make$(cmt_lock_pid)
else
#cmt_local_tagfile_MyIngridLibs = $(INGRID_tag)_MyIngridLibs.make
cmt_local_tagfile_MyIngridLibs = $(bin)$(INGRID_tag)_MyIngridLibs.make
endif

else

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_MyIngridLibs = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile_MyIngridLibs = $(INGRID_tag).make
cmt_local_tagfile_MyIngridLibs = $(bin)$(INGRID_tag).make
endif

endif

-include $(cmt_local_tagfile_MyIngridLibs)

ifdef cmt_MyIngridLibs_has_target_tag

ifdef READONLY
cmt_final_setup_MyIngridLibs = /tmp/CMT_INGRID_MyIngridLibssetup.make
cmt_local_MyIngridLibs_makefile = /tmp/CMT_MyIngridLibs$(cmt_lock_pid).make
else
cmt_final_setup_MyIngridLibs = $(bin)INGRID_MyIngridLibssetup.make
cmt_local_MyIngridLibs_makefile = $(bin)MyIngridLibs.make
endif

else

ifdef READONLY
cmt_final_setup_MyIngridLibs = /tmp/CMT_INGRIDsetup.make
cmt_local_MyIngridLibs_makefile = /tmp/CMT_MyIngridLibs$(cmt_lock_pid).make
else
cmt_final_setup_MyIngridLibs = $(bin)INGRIDsetup.make
cmt_local_MyIngridLibs_makefile = $(bin)MyIngridLibs.make
endif

endif

ifdef READONLY
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
cmt_final_setup = $(bin)INGRIDsetup.make
endif

MyIngridLibs ::


ifdef READONLY
MyIngridLibs ::
	@echo tags=$(tags)
	@echo cmt_local_tagfile=$(cmt_local_tagfile)
endif


dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'MyIngridLibs'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = MyIngridLibs/
MyIngridLibs::
	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

#-- end of make_header ------------------

#-- start of libary_header ---------------

MyIngridLibslibname   = $(bin)$(library_prefix)MyIngridLibs$(library_suffix)
MyIngridLibslib       = $(MyIngridLibslibname).a
MyIngridLibsstamp     = $(bin)MyIngridLibs.stamp
MyIngridLibsshstamp   = $(bin)MyIngridLibs.shstamp

MyIngridLibs :: dirs  MyIngridLibsLIB
	$(echo) "MyIngridLibs ok"

#-- end of libary_header ----------------
#-- start of libary ----------------------

MyIngridLibsLIB :: $(MyIngridLibslib) $(MyIngridLibsshstamp)
	$(echo) "MyIngridLibs : library ok"

$(MyIngridLibslib) :: $(bin)DSTMaker.o $(bin)TINGRID_version.o $(bin)PM_Ch_config.o $(bin)INGRID_BadCh_mapping.o $(bin)ana_MPPC.o $(bin)ana_cosmic.o $(bin)INGRID_Ch_config.o $(bin)INGRID_Dimension.o $(bin)TINGRID_version_Dict.o $(bin)NeutInfoSummary.o $(bin)IngridSimVertexSummary.o $(bin)IngridHitSummary.o $(bin)IngridSimParticleSummaryDict.o $(bin)IngridBasicReconSummaryDict.o $(bin)IngridHitSummaryDict.o $(bin)IngridBasicReconSummary.o $(bin)IngCalib_new.o $(bin)IngridSimVertexSummaryDict.o $(bin)IngridSimHitSummaryDict.o $(bin)IngridSimParticleSummary.o $(bin)IngridTrackSummaryDict.o $(bin)Ingrid1stReducSummaryDict.o $(bin)test_process.o $(bin)test_data.o $(bin)test_read.o $(bin)DumpIngGoodSpill.o $(bin)INGRIDEVENTSUMMARYDict.o $(bin)NeutInfoSummaryDict.o $(bin)IngridSimHitSummary.o $(bin)BeamInfoSummary.o $(bin)test_calib.o $(bin)BeamInfoSummaryDict.o $(bin)tdc.o $(bin)IngridTrackSummary.o $(bin)INGRIDEVENTSUMMARY.o $(bin)IngCalib.o $(bin)IngAddBSD.o $(bin)Ingrid1stReducSummary.o
	$(lib_echo) "static library $@"
	$(lib_silent) $(ar) $(MyIngridLibslib) $?
	$(lib_silent) $(ranlib) $(MyIngridLibslib)
	$(lib_silent) cat /dev/null >$(MyIngridLibsstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

#
# We add one level of dependency upon the true shared library 
# (rather than simply upon the stamp file)
# this is for cases where the shared library has not been built
# while the stamp was created (error??) 
#

$(MyIngridLibslibname).$(shlibsuffix) :: $(MyIngridLibslib) requirements $(use_requirements) $(MyIngridLibsstamps)
	$(lib_echo) "shared library $@"
	$(lib_silent) if test "$(makecmd)"; then QUIET=; else QUIET=1; fi; QUIET=$${QUIET} bin=$(bin) $(make_shlib) "$(tags)" MyIngridLibs $(MyIngridLibs_shlibflags)

$(MyIngridLibsshstamp) :: $(MyIngridLibslibname).$(shlibsuffix)
	$(lib_silent) if test -f $(MyIngridLibslibname).$(shlibsuffix) ; then cat /dev/null >$(MyIngridLibsshstamp) ; fi

MyIngridLibsclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) /bin/rm -f $(bin)DSTMaker.o $(bin)TINGRID_version.o $(bin)PM_Ch_config.o $(bin)INGRID_BadCh_mapping.o $(bin)ana_MPPC.o $(bin)ana_cosmic.o $(bin)INGRID_Ch_config.o $(bin)INGRID_Dimension.o $(bin)TINGRID_version_Dict.o $(bin)NeutInfoSummary.o $(bin)IngridSimVertexSummary.o $(bin)IngridHitSummary.o $(bin)IngridSimParticleSummaryDict.o $(bin)IngridBasicReconSummaryDict.o $(bin)IngridHitSummaryDict.o $(bin)IngridBasicReconSummary.o $(bin)IngCalib_new.o $(bin)IngridSimVertexSummaryDict.o $(bin)IngridSimHitSummaryDict.o $(bin)IngridSimParticleSummary.o $(bin)IngridTrackSummaryDict.o $(bin)Ingrid1stReducSummaryDict.o $(bin)test_process.o $(bin)test_data.o $(bin)test_read.o $(bin)DumpIngGoodSpill.o $(bin)INGRIDEVENTSUMMARYDict.o $(bin)NeutInfoSummaryDict.o $(bin)IngridSimHitSummary.o $(bin)BeamInfoSummary.o $(bin)test_calib.o $(bin)BeamInfoSummaryDict.o $(bin)tdc.o $(bin)IngridTrackSummary.o $(bin)INGRIDEVENTSUMMARY.o $(bin)IngCalib.o $(bin)IngAddBSD.o $(bin)Ingrid1stReducSummary.o
	$(cleanup_silent) cd $(bin); /bin/rm -rf MyIngridLibs_deps MyIngridLibs_dependencies.make

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

install_dir = ${CMTINSTALLAREA}/$(tag)/lib
MyIngridLibsinstallname = $(library_prefix)MyIngridLibs$(library_suffix).$(shlibsuffix)

MyIngridLibs :: MyIngridLibsinstall

install :: MyIngridLibsinstall

MyIngridLibsinstall :: $(install_dir)/$(MyIngridLibsinstallname)
ifdef CMTINSTALLAREA
	$(echo) "installation done"
endif

$(install_dir)/$(MyIngridLibsinstallname) :: $(bin)$(MyIngridLibsinstallname)
ifdef CMTINSTALLAREA
	$(install_silent) $(cmt_install_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(MyIngridLibsinstallname)" \
	    -out "$(install_dir)" \
	    -cmd "$(cmt_installarea_command)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

##MyIngridLibsclean :: MyIngridLibsuninstall

uninstall :: MyIngridLibsuninstall

MyIngridLibsuninstall ::
ifdef CMTINSTALLAREA
	$(cleanup_silent) $(cmt_uninstall_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(MyIngridLibsinstallname)" \
	    -out "$(install_dir)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

#-- end of libary -----------------------
#-- start of dependency ------------------
ifneq ($(MAKECMDGOALS),MyIngridLibsclean)

#$(bin)MyIngridLibs_dependencies.make :: dirs

ifndef QUICK
$(bin)MyIngridLibs_dependencies.make : ../src/DSTMaker.cxx ../src/TINGRID_version.cxx ../src/PM_Ch_config.cxx ../src/INGRID_BadCh_mapping.cxx ../src/ana_MPPC.cxx ../src/ana_cosmic.cxx ../src/INGRID_Ch_config.cxx ../src/INGRID_Dimension.cxx ../dict/TINGRID_version_Dict.cxx ../src/NeutInfoSummary.cc ../src/IngridSimVertexSummary.cc ../src/IngridHitSummary.cc ../src/IngridSimParticleSummaryDict.cc ../src/IngridBasicReconSummaryDict.cc ../src/IngridHitSummaryDict.cc ../src/IngridBasicReconSummary.cc ../src/IngCalib_new.cc ../src/IngridSimVertexSummaryDict.cc ../src/IngridSimHitSummaryDict.cc ../src/IngridSimParticleSummary.cc ../src/IngridTrackSummaryDict.cc ../src/Ingrid1stReducSummaryDict.cc ../src/test_process.cc ../src/test_data.cc ../src/test_read.cc ../src/DumpIngGoodSpill.cc ../src/INGRIDEVENTSUMMARYDict.cc ../src/NeutInfoSummaryDict.cc ../src/IngridSimHitSummary.cc ../src/BeamInfoSummary.cc ../src/test_calib.cc ../src/BeamInfoSummaryDict.cc ../src/tdc.cc ../src/IngridTrackSummary.cc ../src/INGRIDEVENTSUMMARY.cc ../src/IngCalib.cc ../src/IngAddBSD.cc ../src/Ingrid1stReducSummary.cc $(use_requirements) $(cmt_final_setup_MyIngridLibs)
	$(echo) "(MyIngridLibs.make) Rebuilding $@"; \
	  $(build_dependencies) MyIngridLibs -all_sources -out=$@ ../src/DSTMaker.cxx ../src/TINGRID_version.cxx ../src/PM_Ch_config.cxx ../src/INGRID_BadCh_mapping.cxx ../src/ana_MPPC.cxx ../src/ana_cosmic.cxx ../src/INGRID_Ch_config.cxx ../src/INGRID_Dimension.cxx ../dict/TINGRID_version_Dict.cxx ../src/NeutInfoSummary.cc ../src/IngridSimVertexSummary.cc ../src/IngridHitSummary.cc ../src/IngridSimParticleSummaryDict.cc ../src/IngridBasicReconSummaryDict.cc ../src/IngridHitSummaryDict.cc ../src/IngridBasicReconSummary.cc ../src/IngCalib_new.cc ../src/IngridSimVertexSummaryDict.cc ../src/IngridSimHitSummaryDict.cc ../src/IngridSimParticleSummary.cc ../src/IngridTrackSummaryDict.cc ../src/Ingrid1stReducSummaryDict.cc ../src/test_process.cc ../src/test_data.cc ../src/test_read.cc ../src/DumpIngGoodSpill.cc ../src/INGRIDEVENTSUMMARYDict.cc ../src/NeutInfoSummaryDict.cc ../src/IngridSimHitSummary.cc ../src/BeamInfoSummary.cc ../src/test_calib.cc ../src/BeamInfoSummaryDict.cc ../src/tdc.cc ../src/IngridTrackSummary.cc ../src/INGRIDEVENTSUMMARY.cc ../src/IngCalib.cc ../src/IngAddBSD.cc ../src/Ingrid1stReducSummary.cc
endif

#$(MyIngridLibs_dependencies)

-include $(bin)MyIngridLibs_dependencies.make

endif
#-- end of dependency -------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(DSTMaker_cxx_dependencies)

$(bin)$(binobj)DSTMaker.o : $(DSTMaker_cxx_dependencies)
	$(cpp_echo) $(src)DSTMaker.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(DSTMaker_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(DSTMaker_cppflags) $(DSTMaker_cxx_cppflags)  $(src)DSTMaker.cxx

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(TINGRID_version_cxx_dependencies)

$(bin)$(binobj)TINGRID_version.o : $(TINGRID_version_cxx_dependencies)
	$(cpp_echo) $(src)TINGRID_version.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(TINGRID_version_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(TINGRID_version_cppflags) $(TINGRID_version_cxx_cppflags)  $(src)TINGRID_version.cxx

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(PM_Ch_config_cxx_dependencies)

$(bin)$(binobj)PM_Ch_config.o : $(PM_Ch_config_cxx_dependencies)
	$(cpp_echo) $(src)PM_Ch_config.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(PM_Ch_config_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(PM_Ch_config_cppflags) $(PM_Ch_config_cxx_cppflags)  $(src)PM_Ch_config.cxx

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(INGRID_BadCh_mapping_cxx_dependencies)

$(bin)$(binobj)INGRID_BadCh_mapping.o : $(INGRID_BadCh_mapping_cxx_dependencies)
	$(cpp_echo) $(src)INGRID_BadCh_mapping.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(INGRID_BadCh_mapping_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(INGRID_BadCh_mapping_cppflags) $(INGRID_BadCh_mapping_cxx_cppflags)  $(src)INGRID_BadCh_mapping.cxx

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(ana_MPPC_cxx_dependencies)

$(bin)$(binobj)ana_MPPC.o : $(ana_MPPC_cxx_dependencies)
	$(cpp_echo) $(src)ana_MPPC.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(ana_MPPC_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(ana_MPPC_cppflags) $(ana_MPPC_cxx_cppflags)  $(src)ana_MPPC.cxx

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(ana_cosmic_cxx_dependencies)

$(bin)$(binobj)ana_cosmic.o : $(ana_cosmic_cxx_dependencies)
	$(cpp_echo) $(src)ana_cosmic.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(ana_cosmic_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(ana_cosmic_cppflags) $(ana_cosmic_cxx_cppflags)  $(src)ana_cosmic.cxx

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(INGRID_Ch_config_cxx_dependencies)

$(bin)$(binobj)INGRID_Ch_config.o : $(INGRID_Ch_config_cxx_dependencies)
	$(cpp_echo) $(src)INGRID_Ch_config.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(INGRID_Ch_config_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(INGRID_Ch_config_cppflags) $(INGRID_Ch_config_cxx_cppflags)  $(src)INGRID_Ch_config.cxx

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(INGRID_Dimension_cxx_dependencies)

$(bin)$(binobj)INGRID_Dimension.o : $(INGRID_Dimension_cxx_dependencies)
	$(cpp_echo) $(src)INGRID_Dimension.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(INGRID_Dimension_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(INGRID_Dimension_cppflags) $(INGRID_Dimension_cxx_cppflags)  $(src)INGRID_Dimension.cxx

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(TINGRID_version_Dict_cxx_dependencies)

$(bin)$(binobj)TINGRID_version_Dict.o : $(TINGRID_version_Dict_cxx_dependencies)
	$(cpp_echo) ../dict/TINGRID_version_Dict.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(TINGRID_version_Dict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(TINGRID_version_Dict_cppflags) $(TINGRID_version_Dict_cxx_cppflags) -I../dict ../dict/TINGRID_version_Dict.cxx

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(NeutInfoSummary_cc_dependencies)

$(bin)$(binobj)NeutInfoSummary.o : $(NeutInfoSummary_cc_dependencies)
	$(cpp_echo) $(src)NeutInfoSummary.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(NeutInfoSummary_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(NeutInfoSummary_cppflags) $(NeutInfoSummary_cc_cppflags)  $(src)NeutInfoSummary.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridSimVertexSummary_cc_dependencies)

$(bin)$(binobj)IngridSimVertexSummary.o : $(IngridSimVertexSummary_cc_dependencies)
	$(cpp_echo) $(src)IngridSimVertexSummary.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridSimVertexSummary_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridSimVertexSummary_cppflags) $(IngridSimVertexSummary_cc_cppflags)  $(src)IngridSimVertexSummary.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridHitSummary_cc_dependencies)

$(bin)$(binobj)IngridHitSummary.o : $(IngridHitSummary_cc_dependencies)
	$(cpp_echo) $(src)IngridHitSummary.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridHitSummary_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridHitSummary_cppflags) $(IngridHitSummary_cc_cppflags)  $(src)IngridHitSummary.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridSimParticleSummaryDict_cc_dependencies)

$(bin)$(binobj)IngridSimParticleSummaryDict.o : $(IngridSimParticleSummaryDict_cc_dependencies)
	$(cpp_echo) $(src)IngridSimParticleSummaryDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridSimParticleSummaryDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridSimParticleSummaryDict_cppflags) $(IngridSimParticleSummaryDict_cc_cppflags)  $(src)IngridSimParticleSummaryDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridBasicReconSummaryDict_cc_dependencies)

$(bin)$(binobj)IngridBasicReconSummaryDict.o : $(IngridBasicReconSummaryDict_cc_dependencies)
	$(cpp_echo) $(src)IngridBasicReconSummaryDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridBasicReconSummaryDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridBasicReconSummaryDict_cppflags) $(IngridBasicReconSummaryDict_cc_cppflags)  $(src)IngridBasicReconSummaryDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridHitSummaryDict_cc_dependencies)

$(bin)$(binobj)IngridHitSummaryDict.o : $(IngridHitSummaryDict_cc_dependencies)
	$(cpp_echo) $(src)IngridHitSummaryDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridHitSummaryDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridHitSummaryDict_cppflags) $(IngridHitSummaryDict_cc_cppflags)  $(src)IngridHitSummaryDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridBasicReconSummary_cc_dependencies)

$(bin)$(binobj)IngridBasicReconSummary.o : $(IngridBasicReconSummary_cc_dependencies)
	$(cpp_echo) $(src)IngridBasicReconSummary.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridBasicReconSummary_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridBasicReconSummary_cppflags) $(IngridBasicReconSummary_cc_cppflags)  $(src)IngridBasicReconSummary.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngCalib_new_cc_dependencies)

$(bin)$(binobj)IngCalib_new.o : $(IngCalib_new_cc_dependencies)
	$(cpp_echo) $(src)IngCalib_new.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngCalib_new_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngCalib_new_cppflags) $(IngCalib_new_cc_cppflags)  $(src)IngCalib_new.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridSimVertexSummaryDict_cc_dependencies)

$(bin)$(binobj)IngridSimVertexSummaryDict.o : $(IngridSimVertexSummaryDict_cc_dependencies)
	$(cpp_echo) $(src)IngridSimVertexSummaryDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridSimVertexSummaryDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridSimVertexSummaryDict_cppflags) $(IngridSimVertexSummaryDict_cc_cppflags)  $(src)IngridSimVertexSummaryDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridSimHitSummaryDict_cc_dependencies)

$(bin)$(binobj)IngridSimHitSummaryDict.o : $(IngridSimHitSummaryDict_cc_dependencies)
	$(cpp_echo) $(src)IngridSimHitSummaryDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridSimHitSummaryDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridSimHitSummaryDict_cppflags) $(IngridSimHitSummaryDict_cc_cppflags)  $(src)IngridSimHitSummaryDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridSimParticleSummary_cc_dependencies)

$(bin)$(binobj)IngridSimParticleSummary.o : $(IngridSimParticleSummary_cc_dependencies)
	$(cpp_echo) $(src)IngridSimParticleSummary.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridSimParticleSummary_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridSimParticleSummary_cppflags) $(IngridSimParticleSummary_cc_cppflags)  $(src)IngridSimParticleSummary.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridTrackSummaryDict_cc_dependencies)

$(bin)$(binobj)IngridTrackSummaryDict.o : $(IngridTrackSummaryDict_cc_dependencies)
	$(cpp_echo) $(src)IngridTrackSummaryDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridTrackSummaryDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridTrackSummaryDict_cppflags) $(IngridTrackSummaryDict_cc_cppflags)  $(src)IngridTrackSummaryDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(Ingrid1stReducSummaryDict_cc_dependencies)

$(bin)$(binobj)Ingrid1stReducSummaryDict.o : $(Ingrid1stReducSummaryDict_cc_dependencies)
	$(cpp_echo) $(src)Ingrid1stReducSummaryDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(Ingrid1stReducSummaryDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(Ingrid1stReducSummaryDict_cppflags) $(Ingrid1stReducSummaryDict_cc_cppflags)  $(src)Ingrid1stReducSummaryDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(test_process_cc_dependencies)

$(bin)$(binobj)test_process.o : $(test_process_cc_dependencies)
	$(cpp_echo) $(src)test_process.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(test_process_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(test_process_cppflags) $(test_process_cc_cppflags)  $(src)test_process.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(test_data_cc_dependencies)

$(bin)$(binobj)test_data.o : $(test_data_cc_dependencies)
	$(cpp_echo) $(src)test_data.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(test_data_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(test_data_cppflags) $(test_data_cc_cppflags)  $(src)test_data.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(test_read_cc_dependencies)

$(bin)$(binobj)test_read.o : $(test_read_cc_dependencies)
	$(cpp_echo) $(src)test_read.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(test_read_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(test_read_cppflags) $(test_read_cc_cppflags)  $(src)test_read.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(DumpIngGoodSpill_cc_dependencies)

$(bin)$(binobj)DumpIngGoodSpill.o : $(DumpIngGoodSpill_cc_dependencies)
	$(cpp_echo) $(src)DumpIngGoodSpill.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(DumpIngGoodSpill_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(DumpIngGoodSpill_cppflags) $(DumpIngGoodSpill_cc_cppflags)  $(src)DumpIngGoodSpill.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(INGRIDEVENTSUMMARYDict_cc_dependencies)

$(bin)$(binobj)INGRIDEVENTSUMMARYDict.o : $(INGRIDEVENTSUMMARYDict_cc_dependencies)
	$(cpp_echo) $(src)INGRIDEVENTSUMMARYDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(INGRIDEVENTSUMMARYDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(INGRIDEVENTSUMMARYDict_cppflags) $(INGRIDEVENTSUMMARYDict_cc_cppflags)  $(src)INGRIDEVENTSUMMARYDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(NeutInfoSummaryDict_cc_dependencies)

$(bin)$(binobj)NeutInfoSummaryDict.o : $(NeutInfoSummaryDict_cc_dependencies)
	$(cpp_echo) $(src)NeutInfoSummaryDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(NeutInfoSummaryDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(NeutInfoSummaryDict_cppflags) $(NeutInfoSummaryDict_cc_cppflags)  $(src)NeutInfoSummaryDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridSimHitSummary_cc_dependencies)

$(bin)$(binobj)IngridSimHitSummary.o : $(IngridSimHitSummary_cc_dependencies)
	$(cpp_echo) $(src)IngridSimHitSummary.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridSimHitSummary_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridSimHitSummary_cppflags) $(IngridSimHitSummary_cc_cppflags)  $(src)IngridSimHitSummary.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(BeamInfoSummary_cc_dependencies)

$(bin)$(binobj)BeamInfoSummary.o : $(BeamInfoSummary_cc_dependencies)
	$(cpp_echo) $(src)BeamInfoSummary.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(BeamInfoSummary_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(BeamInfoSummary_cppflags) $(BeamInfoSummary_cc_cppflags)  $(src)BeamInfoSummary.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(test_calib_cc_dependencies)

$(bin)$(binobj)test_calib.o : $(test_calib_cc_dependencies)
	$(cpp_echo) $(src)test_calib.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(test_calib_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(test_calib_cppflags) $(test_calib_cc_cppflags)  $(src)test_calib.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(BeamInfoSummaryDict_cc_dependencies)

$(bin)$(binobj)BeamInfoSummaryDict.o : $(BeamInfoSummaryDict_cc_dependencies)
	$(cpp_echo) $(src)BeamInfoSummaryDict.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(BeamInfoSummaryDict_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(BeamInfoSummaryDict_cppflags) $(BeamInfoSummaryDict_cc_cppflags)  $(src)BeamInfoSummaryDict.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(tdc_cc_dependencies)

$(bin)$(binobj)tdc.o : $(tdc_cc_dependencies)
	$(cpp_echo) $(src)tdc.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(tdc_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(tdc_cppflags) $(tdc_cc_cppflags)  $(src)tdc.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngridTrackSummary_cc_dependencies)

$(bin)$(binobj)IngridTrackSummary.o : $(IngridTrackSummary_cc_dependencies)
	$(cpp_echo) $(src)IngridTrackSummary.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngridTrackSummary_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngridTrackSummary_cppflags) $(IngridTrackSummary_cc_cppflags)  $(src)IngridTrackSummary.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(INGRIDEVENTSUMMARY_cc_dependencies)

$(bin)$(binobj)INGRIDEVENTSUMMARY.o : $(INGRIDEVENTSUMMARY_cc_dependencies)
	$(cpp_echo) $(src)INGRIDEVENTSUMMARY.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(INGRIDEVENTSUMMARY_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(INGRIDEVENTSUMMARY_cppflags) $(INGRIDEVENTSUMMARY_cc_cppflags)  $(src)INGRIDEVENTSUMMARY.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngCalib_cc_dependencies)

$(bin)$(binobj)IngCalib.o : $(IngCalib_cc_dependencies)
	$(cpp_echo) $(src)IngCalib.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngCalib_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngCalib_cppflags) $(IngCalib_cc_cppflags)  $(src)IngCalib.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(IngAddBSD_cc_dependencies)

$(bin)$(binobj)IngAddBSD.o : $(IngAddBSD_cc_dependencies)
	$(cpp_echo) $(src)IngAddBSD.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(IngAddBSD_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(IngAddBSD_cppflags) $(IngAddBSD_cc_cppflags)  $(src)IngAddBSD.cc

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

$(bin)MyIngridLibs_dependencies.make : $(Ingrid1stReducSummary_cc_dependencies)

$(bin)$(binobj)Ingrid1stReducSummary.o : $(Ingrid1stReducSummary_cc_dependencies)
	$(cpp_echo) $(src)Ingrid1stReducSummary.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(MyIngridLibs_pp_cppflags) $(lib_MyIngridLibs_pp_cppflags) $(Ingrid1stReducSummary_pp_cppflags) $(use_cppflags) $(MyIngridLibs_cppflags) $(lib_MyIngridLibs_cppflags) $(Ingrid1stReducSummary_cppflags) $(Ingrid1stReducSummary_cc_cppflags)  $(src)Ingrid1stReducSummary.cc

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: MyIngridLibsclean
	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(MyIngridLibs.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
else
.DEFAULT::
	$(echo) "(MyIngridLibs.make) PEDANTIC: $@: No rule for such target" >&2
	if test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_MyIngridLibs)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(MyIngridLibs.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr index $@ '/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(MyIngridLibs.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(MyIngridLibs.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

MyIngridLibsclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(binobj)$(library_prefix)MyIngridLibs$(library_suffix).a $(binobj)$(library_prefix)MyIngridLibs$(library_suffix).s? $(binobj)MyIngridLibs.stamp $(binobj)MyIngridLibs.shstamp
#-- end of cleanup_library ---------------

