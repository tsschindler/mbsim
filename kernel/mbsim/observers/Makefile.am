noinst_LTLIBRARIES = libobserver.la
libobserver_la_SOURCES = kinematics_observer.cc \
			 coordinates_observer.cc \
			 frame_observer.cc \
			 rigid_body_observer.cc \
			 rigid_body_group_observer.cc

observerincludedir = $(includedir)/mbsim/observers

libobserver_la_LIBADD = $(DEPS_LIBS) $(OPENMBVCPPINTERFACE_LIBS)
libobserver_la_CPPFLAGS = -I$(top_srcdir) $(DEPS_CFLAGS) $(OPENMBVCPPINTERFACE_CFLAGS)

observerinclude_HEADERS = kinematics_observer.h \
			  coordinates_observer.h \
			  frame_observer.h \
			  rigid_body_observer.h \
			  rigid_body_group_observer.h
