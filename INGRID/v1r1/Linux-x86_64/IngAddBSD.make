#-- start of make_header -----------------

#====================================
#  Application IngAddBSD
#
#   Generated Sun Mar  3 01:44:44 2013  by bquilain
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_IngAddBSD_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_IngAddBSD_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_IngAddBSD

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_IngAddBSD = /tmp/CMT_$(INGRID_tag)_IngAddBSD.make$(cmt_lock_pid)
else
#cmt_local_tagfile_IngAddBSD = $(INGRID_tag)_IngAddBSD.make
cmt_local_tagfile_IngAddBSD = $(bin)$(INGRID_tag)_IngAddBSD.make
endif

else

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_IngAddBSD = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile_IngAddBSD = $(INGRID_tag).make
cmt_local_tagfile_IngAddBSD = $(bin)$(INGRID_tag).make
endif

endif

-include $(cmt_local_tagfile_IngAddBSD)

ifdef cmt_IngAddBSD_has_target_tag

ifdef READONLY
cmt_final_setup_IngAddBSD = /tmp/CMT_INGRID_IngAddBSDsetup.make
cmt_local_IngAddBSD_makefile = /tmp/CMT_IngAddBSD$(cmt_lock_pid).make
else
cmt_final_setup_IngAddBSD = $(bin)INGRID_IngAddBSDsetup.make
cmt_local_IngAddBSD_makefile = $(bin)IngAddBSD.make
endif

else

ifdef READONLY
cmt_final_setup_IngAddBSD = /tmp/CMT_INGRIDsetup.make
cmt_local_IngAddBSD_makefile = /tmp/CMT_IngAddBSD$(cmt_lock_pid).make
else
cmt_final_setup_IngAddBSD = $(bin)INGRIDsetup.make
cmt_local_IngAddBSD_makefile = $(bin)IngAddBSD.make
endif

endif

ifdef READONLY
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
cmt_final_setup = $(bin)INGRIDsetup.make
endif

IngAddBSD ::


ifdef READONLY
IngAddBSD ::
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
	$(echo) 'IngAddBSD'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = IngAddBSD/
IngAddBSD::
	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

#-- end of make_header ------------------

#-- start of application_header

IngAddBSD :: dirs  $(bin)IngAddBSD${application_suffix}
	$(echo) "IngAddBSD ok"

#-- end of application_header
#-- start of application

$(bin)IngAddBSD${application_suffix} :: $(bin)IngAddBSD.o $(use_stamps) $(IngAddBSDstamps) requirements $(use_requirements)
	$(link_echo) "application $@"
	$(link_silent) $(cpplink) -o $(@).new $(bin)IngAddBSD.o $(cmt_installarea_linkopts) $(IngAddBSD_use_linkopts) $(IngAddBSDlinkopts) && mv -f $(@).new $(@)

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

install_dir = ${CMTINSTALLAREA}/$(tag)/bin
IngAddBSDinstallname = IngAddBSD${application_suffix}

IngAddBSD :: IngAddBSDinstall

install :: IngAddBSDinstall

IngAddBSDinstall :: $(install_dir)/$(IngAddBSDinstallname)
ifdef CMTINSTALLAREA
	$(echo) "installation done"
endif

$(install_dir)/$(IngAddBSDinstallname) :: $(bin)$(IngAddBSDinstallname)
ifdef CMTINSTALLAREA
	$(install_silent) $(cmt_install_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(IngAddBSDinstallname)" \
	    -out "$(install_dir)" \
	    -cmd "$(cmt_installarea_command)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

##IngAddBSDclean :: IngAddBSDuninstall

uninstall :: IngAddBSDuninstall

IngAddBSDuninstall ::
ifdef CMTINSTALLAREA
	$(cleanup_silent) $(cmt_uninstall_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(IngAddBSDinstallname)" \
	    -out "$(install_dir)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

#	@echo "------> (IngAddBSD.make) Removing installed files"
#-- end of application
#-- start of dependency ------------------
ifneq ($(MAKECMDGOALS),IngAddBSDclean)

#$(bin)IngAddBSD_dependencies.make :: dirs

ifndef QUICK
$(bin)IngAddBSD_dependencies.make : ../app/IngAddBSD.cc $(use_requirements) $(cmt_final_setup_IngAddBSD)
	$(echo) "(IngAddBSD.make) Rebuilding $@"; \
	  $(build_dependencies) IngAddBSD -all_sources -out=$@ ../app/IngAddBSD.cc
endif

#$(IngAddBSD_dependencies)

-include $(bin)IngAddBSD_dependencies.make

endif
#-- end of dependency -------------------
#-- start of cpp ------

$(bin)IngAddBSD_dependencies.make : $(IngAddBSD_cc_dependencies)

$(bin)$(binobj)IngAddBSD.o : $(IngAddBSD_cc_dependencies)
	$(cpp_echo) ../app/IngAddBSD.cc
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(IngAddBSD_pp_cppflags) $(app_IngAddBSD_pp_cppflags) $(IngAddBSD_pp_cppflags) $(use_cppflags) $(IngAddBSD_cppflags) $(app_IngAddBSD_cppflags) $(IngAddBSD_cppflags) $(IngAddBSD_cc_cppflags) -I../app ../app/IngAddBSD.cc

#-- end of cpp ------
#-- start of cleanup_header --------------

clean :: IngAddBSDclean
	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(IngAddBSD.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
else
.DEFAULT::
	$(echo) "(IngAddBSD.make) PEDANTIC: $@: No rule for such target" >&2
	if test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_IngAddBSD)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(IngAddBSD.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr index $@ '/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(IngAddBSD.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(IngAddBSD.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

IngAddBSDclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_application ------
	$(cleanup_echo) IngAddBSD${application_suffix}
	-$(cleanup_silent) cd $(bin); /bin/rm -f IngAddBSD${application_suffix}

#	@echo "------> (IngAddBSD.make) Removing application files"
#-- end of cleanup_application ------
#-- start of cleanup_objects ------
	$(cleanup_echo) objects
	-$(cleanup_silent) /bin/rm -f $(bin)IngAddBSD.o
	-$(cleanup_silent) cd $(bin); /bin/rm -rf IngAddBSD_deps IngAddBSD_dependencies.make
#-- end of cleanup_objects ------
