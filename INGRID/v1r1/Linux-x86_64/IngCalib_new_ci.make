#-- start of make_header -----------------

#====================================
#  Application IngCalib_new_ci
#
#   Generated Tue Mar 19 01:14:14 2013  by bquilain
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_IngCalib_new_ci_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_IngCalib_new_ci_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_IngCalib_new_ci

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_IngCalib_new_ci = /tmp/CMT_$(INGRID_tag)_IngCalib_new_ci.make$(cmt_lock_pid)
else
#cmt_local_tagfile_IngCalib_new_ci = $(INGRID_tag)_IngCalib_new_ci.make
cmt_local_tagfile_IngCalib_new_ci = $(bin)$(INGRID_tag)_IngCalib_new_ci.make
endif

else

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_IngCalib_new_ci = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile_IngCalib_new_ci = $(INGRID_tag).make
cmt_local_tagfile_IngCalib_new_ci = $(bin)$(INGRID_tag).make
endif

endif

-include $(cmt_local_tagfile_IngCalib_new_ci)

ifdef cmt_IngCalib_new_ci_has_target_tag

ifdef READONLY
cmt_final_setup_IngCalib_new_ci = /tmp/CMT_INGRID_IngCalib_new_cisetup.make
cmt_local_IngCalib_new_ci_makefile = /tmp/CMT_IngCalib_new_ci$(cmt_lock_pid).make
else
cmt_final_setup_IngCalib_new_ci = $(bin)INGRID_IngCalib_new_cisetup.make
cmt_local_IngCalib_new_ci_makefile = $(bin)IngCalib_new_ci.make
endif

else

ifdef READONLY
cmt_final_setup_IngCalib_new_ci = /tmp/CMT_INGRIDsetup.make
cmt_local_IngCalib_new_ci_makefile = /tmp/CMT_IngCalib_new_ci$(cmt_lock_pid).make
else
cmt_final_setup_IngCalib_new_ci = $(bin)INGRIDsetup.make
cmt_local_IngCalib_new_ci_makefile = $(bin)IngCalib_new_ci.make
endif

endif

ifdef READONLY
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
cmt_final_setup = $(bin)INGRIDsetup.make
endif

IngCalib_new_ci ::


ifdef READONLY
IngCalib_new_ci ::
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
	$(echo) 'IngCalib_new_ci'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = IngCalib_new_ci/
IngCalib_new_ci::
	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

#-- end of make_header ------------------

#-- start of application_header

IngCalib_new_ci :: dirs  $(bin)IngCalib_new_ci${application_suffix}
	$(echo) "IngCalib_new_ci ok"

#-- end of application_header
#-- start of application

$(bin)IngCalib_new_ci${application_suffix} :: $(bin)IngCalib_new_ci.o $(use_stamps) $(IngCalib_new_cistamps) requirements $(use_requirements)
	$(link_echo) "application $@"
	$(link_silent) $(cpplink) -o $(@).new $(bin)IngCalib_new_ci.o $(cmt_installarea_linkopts) $(IngCalib_new_ci_use_linkopts) $(IngCalib_new_cilinkopts) && mv -f $(@).new $(@)

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

install_dir = ${CMTINSTALLAREA}/$(tag)/bin
IngCalib_new_ciinstallname = IngCalib_new_ci${application_suffix}

IngCalib_new_ci :: IngCalib_new_ciinstall

install :: IngCalib_new_ciinstall

IngCalib_new_ciinstall :: $(install_dir)/$(IngCalib_new_ciinstallname)
ifdef CMTINSTALLAREA
	$(echo) "installation done"
endif

$(install_dir)/$(IngCalib_new_ciinstallname) :: $(bin)$(IngCalib_new_ciinstallname)
ifdef CMTINSTALLAREA
	$(install_silent) $(cmt_install_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(IngCalib_new_ciinstallname)" \
	    -out "$(install_dir)" \
	    -cmd "$(cmt_installarea_command)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

##IngCalib_new_ciclean :: IngCalib_new_ciuninstall

uninstall :: IngCalib_new_ciuninstall

IngCalib_new_ciuninstall ::
ifdef CMTINSTALLAREA
	$(cleanup_silent) $(cmt_uninstall_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(IngCalib_new_ciinstallname)" \
	    -out "$(install_dir)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

#	@echo "------> (IngCalib_new_ci.make) Removing installed files"
#-- end of application
#-- start of dependency ------------------
ifneq ($(MAKECMDGOALS),IngCalib_new_ciclean)

#$(bin)IngCalib_new_ci_dependencies.make :: dirs

ifndef QUICK
$(bin)IngCalib_new_ci_dependencies.make : ../app/IngCalib_new_ci.cc $(use_requirements) $(cmt_final_setup_IngCalib_new_ci)
	$(echo) "(IngCalib_new_ci.make) Rebuilding $@"; \
	  $(build_dependencies) IngCalib_new_ci -all_sources -out=$@ ../app/IngCalib_new_ci.cc
endif

#$(IngCalib_new_ci_dependencies)

-include $(bin)IngCalib_new_ci_dependencies.make

endif
#-- end of dependency -------------------
#-- start of cpp ------

$(bin)IngCalib_new_ci_dependencies.make : $(IngCalib_new_ci_cc_dependencies)

$(bin)$(binobj)IngCalib_new_ci.o : $(IngCalib_new_ci_cc_dependencies)
	$(cpp_echo) ../app/IngCalib_new_ci.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(IngCalib_new_ci_pp_cppflags) $(app_IngCalib_new_ci_pp_cppflags) $(IngCalib_new_ci_pp_cppflags) $(use_cppflags) $(IngCalib_new_ci_cppflags) $(app_IngCalib_new_ci_cppflags) $(IngCalib_new_ci_cppflags) $(IngCalib_new_ci_cc_cppflags) -I../app ../app/IngCalib_new_ci.cc

#-- end of cpp ------
#-- start of cleanup_header --------------

clean :: IngCalib_new_ciclean
	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(IngCalib_new_ci.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
else
.DEFAULT::
	$(echo) "(IngCalib_new_ci.make) PEDANTIC: $@: No rule for such target" >&2
	if test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_IngCalib_new_ci)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(IngCalib_new_ci.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr index $@ '/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(IngCalib_new_ci.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(IngCalib_new_ci.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

IngCalib_new_ciclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_application ------
	$(cleanup_echo) IngCalib_new_ci${application_suffix}
	-$(cleanup_silent) cd $(bin); /bin/rm -f IngCalib_new_ci${application_suffix}

#	@echo "------> (IngCalib_new_ci.make) Removing application files"
#-- end of cleanup_application ------
#-- start of cleanup_objects ------
	$(cleanup_echo) objects
	-$(cleanup_silent) /bin/rm -f $(bin)IngCalib_new_ci.o
	-$(cleanup_silent) cd $(bin); /bin/rm -rf IngCalib_new_ci_deps IngCalib_new_ci_dependencies.make
#-- end of cleanup_objects ------
