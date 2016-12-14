#-- start of make_header -----------------

#====================================
#  Application Calc_MPPC_new_ci
#
#   Generated Sun Aug 17 11:05:41 2014  by bquilain
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_Calc_MPPC_new_ci_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_Calc_MPPC_new_ci_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_Calc_MPPC_new_ci

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_Calc_MPPC_new_ci = /tmp/CMT_$(INGRID_tag)_Calc_MPPC_new_ci.make$(cmt_lock_pid)
else
#cmt_local_tagfile_Calc_MPPC_new_ci = $(INGRID_tag)_Calc_MPPC_new_ci.make
cmt_local_tagfile_Calc_MPPC_new_ci = $(bin)$(INGRID_tag)_Calc_MPPC_new_ci.make
endif

else

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_Calc_MPPC_new_ci = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile_Calc_MPPC_new_ci = $(INGRID_tag).make
cmt_local_tagfile_Calc_MPPC_new_ci = $(bin)$(INGRID_tag).make
endif

endif

-include $(cmt_local_tagfile_Calc_MPPC_new_ci)

ifdef cmt_Calc_MPPC_new_ci_has_target_tag

ifdef READONLY
cmt_final_setup_Calc_MPPC_new_ci = /tmp/CMT_INGRID_Calc_MPPC_new_cisetup.make
cmt_local_Calc_MPPC_new_ci_makefile = /tmp/CMT_Calc_MPPC_new_ci$(cmt_lock_pid).make
else
cmt_final_setup_Calc_MPPC_new_ci = $(bin)INGRID_Calc_MPPC_new_cisetup.make
cmt_local_Calc_MPPC_new_ci_makefile = $(bin)Calc_MPPC_new_ci.make
endif

else

ifdef READONLY
cmt_final_setup_Calc_MPPC_new_ci = /tmp/CMT_INGRIDsetup.make
cmt_local_Calc_MPPC_new_ci_makefile = /tmp/CMT_Calc_MPPC_new_ci$(cmt_lock_pid).make
else
cmt_final_setup_Calc_MPPC_new_ci = $(bin)INGRIDsetup.make
cmt_local_Calc_MPPC_new_ci_makefile = $(bin)Calc_MPPC_new_ci.make
endif

endif

ifdef READONLY
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
cmt_final_setup = $(bin)INGRIDsetup.make
endif

Calc_MPPC_new_ci ::


ifdef READONLY
Calc_MPPC_new_ci ::
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
	$(echo) 'Calc_MPPC_new_ci'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = Calc_MPPC_new_ci/
Calc_MPPC_new_ci::
	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

#-- end of make_header ------------------

#-- start of application_header

Calc_MPPC_new_ci :: dirs  $(bin)Calc_MPPC_new_ci${application_suffix}
	$(echo) "Calc_MPPC_new_ci ok"

#-- end of application_header
#-- start of application

$(bin)Calc_MPPC_new_ci${application_suffix} :: $(bin)Calc_MPPC_new_ci.o $(use_stamps) $(Calc_MPPC_new_cistamps) requirements $(use_requirements)
	$(link_echo) "application $@"
	$(link_silent) $(cpplink) -o $(@).new $(bin)Calc_MPPC_new_ci.o $(cmt_installarea_linkopts) $(Calc_MPPC_new_ci_use_linkopts) $(Calc_MPPC_new_cilinkopts) && mv -f $(@).new $(@)

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

install_dir = ${CMTINSTALLAREA}/$(tag)/bin
Calc_MPPC_new_ciinstallname = Calc_MPPC_new_ci${application_suffix}

Calc_MPPC_new_ci :: Calc_MPPC_new_ciinstall

install :: Calc_MPPC_new_ciinstall

Calc_MPPC_new_ciinstall :: $(install_dir)/$(Calc_MPPC_new_ciinstallname)
ifdef CMTINSTALLAREA
	$(echo) "installation done"
endif

$(install_dir)/$(Calc_MPPC_new_ciinstallname) :: $(bin)$(Calc_MPPC_new_ciinstallname)
ifdef CMTINSTALLAREA
	$(install_silent) $(cmt_install_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(Calc_MPPC_new_ciinstallname)" \
	    -out "$(install_dir)" \
	    -cmd "$(cmt_installarea_command)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

##Calc_MPPC_new_ciclean :: Calc_MPPC_new_ciuninstall

uninstall :: Calc_MPPC_new_ciuninstall

Calc_MPPC_new_ciuninstall ::
ifdef CMTINSTALLAREA
	$(cleanup_silent) $(cmt_uninstall_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(Calc_MPPC_new_ciinstallname)" \
	    -out "$(install_dir)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

#	@echo "------> (Calc_MPPC_new_ci.make) Removing installed files"
#-- end of application
#-- start of dependency ------------------
ifneq ($(MAKECMDGOALS),Calc_MPPC_new_ciclean)

#$(bin)Calc_MPPC_new_ci_dependencies.make :: dirs

ifndef QUICK
$(bin)Calc_MPPC_new_ci_dependencies.make : ../app/Calc_MPPC_new_ci.cxx $(use_requirements) $(cmt_final_setup_Calc_MPPC_new_ci)
	$(echo) "(Calc_MPPC_new_ci.make) Rebuilding $@"; \
	  $(build_dependencies) Calc_MPPC_new_ci -all_sources -out=$@ ../app/Calc_MPPC_new_ci.cxx
endif

#$(Calc_MPPC_new_ci_dependencies)

-include $(bin)Calc_MPPC_new_ci_dependencies.make

endif
#-- end of dependency -------------------
#-- start of cpp ------

$(bin)Calc_MPPC_new_ci_dependencies.make : $(Calc_MPPC_new_ci_cxx_dependencies)

$(bin)$(binobj)Calc_MPPC_new_ci.o : $(Calc_MPPC_new_ci_cxx_dependencies)
	$(cpp_echo) ../app/Calc_MPPC_new_ci.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(Calc_MPPC_new_ci_pp_cppflags) $(app_Calc_MPPC_new_ci_pp_cppflags) $(Calc_MPPC_new_ci_pp_cppflags) $(use_cppflags) $(Calc_MPPC_new_ci_cppflags) $(app_Calc_MPPC_new_ci_cppflags) $(Calc_MPPC_new_ci_cppflags) $(Calc_MPPC_new_ci_cxx_cppflags) -I../app ../app/Calc_MPPC_new_ci.cxx

#-- end of cpp ------
#-- start of cleanup_header --------------

clean :: Calc_MPPC_new_ciclean
	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(Calc_MPPC_new_ci.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
else
.DEFAULT::
	$(echo) "(Calc_MPPC_new_ci.make) PEDANTIC: $@: No rule for such target" >&2
	if test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_Calc_MPPC_new_ci)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(Calc_MPPC_new_ci.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr index $@ '/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(Calc_MPPC_new_ci.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(Calc_MPPC_new_ci.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

Calc_MPPC_new_ciclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_application ------
	$(cleanup_echo) Calc_MPPC_new_ci${application_suffix}
	-$(cleanup_silent) cd $(bin); /bin/rm -f Calc_MPPC_new_ci${application_suffix}

#	@echo "------> (Calc_MPPC_new_ci.make) Removing application files"
#-- end of cleanup_application ------
#-- start of cleanup_objects ------
	$(cleanup_echo) objects
	-$(cleanup_silent) /bin/rm -f $(bin)Calc_MPPC_new_ci.o
	-$(cleanup_silent) cd $(bin); /bin/rm -rf Calc_MPPC_new_ci_deps Calc_MPPC_new_ci_dependencies.make
#-- end of cleanup_objects ------
