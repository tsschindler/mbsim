noinst_LTLIBRARIES = libintegrators.la
libintegrators_la_SOURCES = integrator.cc\
					integratorDAE.cc\
			    lsodar_integrator.cc \
			    lsoder_integrator.cc \
			    lsode_integrator.cc \
			    rksuite_integrator.cc \
			    dopri5_integrator.cc \
			    dop853_integrator.cc \
			    odex_integrator.cc \
			    quasi_static_integrator.cc \
			    radau5_integrator.cc \
			    time_stepping_integrator.cc \
			    time_stepping_ssc_integrator.cc \
			    auto_time_stepping_ssc_integrator.cc \
			    hets2_integrator.cc \
			    theta_time_stepping_integrator.cc \
			    euler_explicit_integrator.cc \
			    daskr_integrator.cc\
			    phem56_integrator.cc\
			    radau5dae_integrator.cc \
			    generalized_alpha_integrator.cc\
			    fortran/opkda1.f\
			    fortran/opkda2.f\
			    fortran/opkdmain.f\
			    fortran/dopri5.f\
			    fortran/dop853.f\
			    fortran/odex.f\
			    fortran/radau5.f\
			    fortran/dc_lapack.f\
			    fortran/rksuite.f\
			    fortran/daux.f\
			    fortran/decsol.f\
			    fortran/linsp.f\
			    fortran/phem56.f\
			    fortran/DASKR/preconds/dbanpre.f\
			    fortran/DASKR/preconds/dilupre.f\
			    fortran/DASKR/preconds/drbdpre.f\
			    fortran/DASKR/preconds/drbgpre.f\
			    fortran/DASKR/preconds/dsparsk.f\
			    fortran/DASKR/solver/ddaskr.f
			    


#			    daspk_integrator.cc			    			    
#			    fortran/DASKR/solver/dlinpk.f
			     
#daspk_integrator.cc \
#fortran/ddaspk.f\
#fortran/daux.f

#fortran/mexax.f
#mexax_integrator.cc 

libintegrators_la_CPPFLAGS = -I$(top_srcdir) $(DEPS_CFLAGS) $(OPENMBVCPPINTERFACE_CFLAGS)
libintegrators_la_LIBADD = $(DEPS_LIBS) $(OPENMBVCPPINTERFACE_LIBS)

integratorsincludedir = $(includedir)/mbsim/integrators

integratorsinclude_HEADERS = integrator.h \
					 integratorDAE.h\
			     integrators.h \
			     time_stepping_integrator.h \
			     time_stepping_ssc_integrator.h \
			     auto_time_stepping_ssc_integrator.h \
			     hets2_integrator.h \
			     theta_time_stepping_integrator.h \
			     euler_explicit_integrator.h \
			     lsodar_integrator.h \
			     lsoder_integrator.h \
			     lsode_integrator.h \
			     rksuite_integrator.h \
			     dopri5_integrator.h \
			     dop853_integrator.h \
			     odex_integrator.h \
			     quasi_static_integrator.h \
			     radau5_integrator.h \
			     daskr_integrator.h \
			     phem56_integrator.h \
			     radau5dae_integrator.h \
			     generalized_alpha_integrator.h\
			     fortran/fortran_wrapper.h
#daspk_integrator.h

# mexax_integrator.h 

