#-- start of make_header -----------------

#====================================
#  Application IngCalib_woCalib
#
#   Generated Fri May 10 19:28:36 2013  by bquilain
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_IngCalib_woCalib_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_IngCalib_woCalib_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_IngCalib_woCalib

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_IngCalib_woCalib = /tmp/CMT_$(INGRID_tag)_IngCalib_woCalib.make$(cmt_lock_pid)
else
#cmt_local_tagfile_IngCalib_woCalib = $(INGRID_tag)_IngCalib_woCalib.make
cmt_local_tagfile_IngCalib_woCalib = $(bin)$(INGRID_tag)_IngCalib_woCalib.make
endif

else

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_IngCalib_woCalib = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile_IngCalib_woCalib = $(INGRID_tag).make
cmt_local_tagfile_IngCalib_woCalib = $(bin)$(INGRID_tag).make
endif

endif

-include $(cmt_local_tagfile_IngCalib_woCalib)

ifdef cmt_IngCalib_woCalib_has_target_tag

ifdef READONLY
cmt_final_setup_IngCalib_woCalib = /tmp/CMT_INGRID_IngCalib_woCalibsetup.make
cmt_local_IngCalib_woCalib_makefile = /tmp/CMT_IngCalib_woCalib$(cmt_lock_pid).make
else
cmt_final_setup_IngCalib_woCalib = $(bin)INGRID_IngCalib_woCalibsetup.make
cmt_local_IngCalib_woCalib_makefile = $(bin)IngCalib_woCalib.make
endif

else

ifdef READONLY
cmt_final_setup_IngCalib_woCalib = /tmp/CMT_INGRIDsetup.make
cmt_local_IngCalib_woCalib_makefile = /tmp/CMT_IngCalib_woCalib$(cmt_lock_pid).make
else
cmt_final_setup_IngCalib_woCalib = $(bin)INGRIDsetup.make
cmt_local_IngCalib_woCalib_makefile = $(bin)IngCalib_woCalib.make
endif

endif

ifdef READONLY
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
cmt_final_setup = $(bin)INGRIDsetup.make
endif

IngCalib_woCalib ::


ifdef READONLY
IngCalib_woCalib ::
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
	$(echo) 'IngCalib_woCalib'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = IngCalib_woCalib/
IngCalib_woCalib::
	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

#-- end of make_header ------------------

#-- start of application_header

IngCalib_woCalib :: dirs  $(bin)IngCalib_woCalib${application_suffix}
	$(echo) "IngCalib_woCalib ok"

#-- end of application_header
#-- start of application

$(bin)IngCalib_woCalib${application_suffix} :: $(bin)IngCalib_woCalib.o $(use_stamps) $(IngCalib_woCalibstamps) requirements $(use_requirements)
	$(link_echo) "application $@"
	$(link_silent) $(cpplink) -o $(@).new $(bin)IngCalib_woCalib.o $(cmt_installarea_linkopts) $(IngCalib_woCalib_use_linkopts) $(IngCalib_woCaliblinkopts) && mv -f $(@).new $(@)

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

install_dir = ${CMTINSTALLAREA}/$(tag)/bin
IngCalib_woCalibinstallname = IngCalib_woCalib${application_suffix}

IngCalib_woCalib :: IngCalib_woCalibinstall

install :: IngCalib_woCalibinstall

IngCalib_woCalibinstall :: $(install_dir)/$(IngCalib_woCalibinstallname)
ifdef CMTINSTALLAREA
	$(echo) "installation done"
endif

$(install_dir)/$(IngCalib_woCalibinstallname) :: $(bin)$(IngCalib_woCalibinstallname)
ifdef CMTINSTALLAREA
	$(install_silent) $(cmt_install_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(IngCalib_woCalibinstallname)" \
	    -out "$(install_dir)" \
	    -cmd "$(cmt_installarea_command)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

##IngCalib_woCalibclean :: IngCalib_woCalibuninstall

uninstall :: IngCalib_woCalibuninstall

IngCalib_woCalibuninstall ::
ifdef CMTINSTALLAREA
	$(cleanup_silent) $(cmt_uninstall_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(IngCalib_woCalibinstallname)" \
	    -out "$(install_dir)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

#	@echo "------> (IngCalib_woCalib.make) Removing installed files"
#-- end of application
#-- start of dependency ------------------
ifneq ($(MAKECMDGOALS),IngCalib_woCalibclean)

#$(bin)IngCalib_woCalib_dependencies.make :: dirs

ifndef QUICK
$(bin)IngCalib_woCalib_dependencies.make : ../app/IngCalib_woCalib.cc $(use_requirements) $(cmt_final_setup_IngCalib_woCalib)
	$(echo) "(IngCalib_woCalib.make) Rebuilding $@"; \
	  $(build_dependencies) IngCalib_woCalib -all_sources -out=$@ ../app/IngCalib_woCalib.cc
endif

#$(IngCalib_woCalib_dependencies)

-include $(bin)IngCalib_woCalib_dependencies.make

endif
#-- end of dependency -------------------
#-- start of cpp ------

$(bin)IngCalib_woCalib_dependencies.make : $(IngCalib_woCalib_cc_dependencies)

$(bin)$(binobj)IngCalib_woCalib.o : $(IngCalib_woCalib_cc_dependencies)
	$(cpp_echo) ../app/IngCalib_woCalib.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(IngCalib_woCalib_pp_cppflags) $(app_IngCalib_woCalib_pp_cppflags) $(IngCalib_woCalib_pp_cppflags) $(use_cppflags) $(IngCalib_woCalib_cppflags) $(app_IngCalib_woCalib_cppflags) $(IngCalib_woCalib_cppflags) $(IngCalib_woCalib_cc_cppflags) -I../app ../app/IngCalib_woCalib.cc

#-- end of cpp ------
#-- start of cleanup_header --------------

clean :: IngCalib_woCalibclean
	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(IngCalib_woCalib.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
else
.DEFAULT::
	$(echo) "(IngCalib_woCalib.make) PEDANTIC: $@: No rule for such target" >&2
	if test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_IngCalib_woCalib)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(IngCalib_woCalib.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr index $@ '/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(IngCalib_woCalib.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(IngCalib_woCalib.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

IngCalib_woCalibclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_application ------
	$(cleanup_echo) IngCalib_woCalib${application_suffix}
	-$(cleanup_silent) cd $(bin); /bin/rm -f IngCalib_woCalib${application_suffix}

#	@echo "------> (IngCalib_woCalib.make) Removing application files"
#-- end of cleanup_application ------
#-- start of cleanup_objects ------
	$(cleanup_echo) objects
	-$(cleanup_silent) /bin/rm -f $(bin)IngCalib_woCalib.o
	-$(cleanup_silent) cd $(bin); /bin/rm -rf IngCalib_woCalib_deps IngCalib_woCalib_dependencies.make
#-- end of cleanup_objects ------
