if HAVE_FORTRAN
AM_CPPFLAGS = -I$(top_builddir) -I$(top_builddir)/include -I$(top_srcdir)/include
LDADD = libarprec_f_main.a libarprecmod.a $(top_builddir)/src/libarprec.a $(FCLIBS)

lib_LIBRARIES = libarprecmod.a libarprec_f_main.a
libarprecmod_a_SOURCES = f_mp.cpp mp_mod.f mp_modm.f mp_modx.f
libarprec_f_main_a_SOURCES = main.cpp

if UPCASE_MODULE
MPCOMPLEXMOD = MPCOMPLEXMOD.$(module_ext)
MPDEFMOD = MPDEFMOD.$(module_ext)
MPGENMOD = MPGENMOD.$(module_ext)
MPINTMOD = MPINTMOD.$(module_ext)
MPMODULE = MPMODULE.$(module_ext)
MPMODULEM = MPMODULEM.$(module_ext)
MPMODULEX = MPMODULEX.$(module_ext)
MPREALMOD = MPREALMOD.$(module_ext)
else
MPCOMPLEXMOD = mpcomplexmod.$(module_ext)
MPDEFMOD = mpdefmod.$(module_ext)
MPGENMOD = mpgenmod.$(module_ext)
MPINTMOD = mpintmod.$(module_ext)
MPMODULE = mpmodule.$(module_ext)
MPMODULEM = mpmodulem.$(module_ext)
MPMODULEX = mpmodulex.$(module_ext)
MPREALMOD = mprealmod.$(module_ext)
endif

MODULES = $(MPCOMPLEXMOD) $(MPDEFMOD) $(MPGENMOD) $(MPINTMOD) \
	  $(MPMODULE) $(MPMODULEM) $(MPMODULEX) $(MPREALMOD)

$(MPDEFMOD) $(MPINTMOD) $(MPREALMOD) $(MPCOMPLEXMOD) $(MPGENMOD) $(MPMODULE): mp_mod.o

mp_modm.o: mp_mod.o
$(MPMODULEM): mp_modm.o

mp_modx.o: mp_mod.o
$(MPMODULEX): mp_modx.o

pkglib_DATA = $(MODULES)

DEMO = pslq1$(EXEEXT) pslq2$(EXEEXT) pslq3$(EXEEXT) pslqm1$(EXEEXT) \
       pslqm2$(EXEEXT) pslqm3$(EXEEXT) quadgs$(EXEEXT) quaderf$(EXEEXT) \
       quadts$(EXEEXT) roots$(EXEEXT)
EXTRA_PROGRAMS = $(DEMO)
EXTRA_DIST = Makefile.sample

pslq1_SOURCES = tpslq1.f
nodist_pslq1_SOURCES = second.f
pslq1_LINK = $(CXXLINK)

pslq2_SOURCES = tpslq2.f
nodist_pslq2_SOURCES = second.f
pslq2_LINK = $(CXXLINK)

pslq3_SOURCES = tpslq3.f
nodist_pslq3_SOURCES = second.f
pslq3_LINK = $(CXXLINK)

pslqm1_SOURCES = tpslqm1.f
nodist_pslqm1_SOURCES = second.f
pslqm1_LINK = $(CXXLINK)

pslqm2_SOURCES = tpslqm2.f
nodist_pslqm2_SOURCES = second.f
pslqm2_LINK = $(CXXLINK)

pslqm3_SOURCES = tpslqm3.f
nodist_pslqm3_SOURCES = second.f
pslqm3_LINK = $(CXXLINK)

quadgs_SOURCES = tquadgs.f
nodist_quadgs_SOURCES = second.f
quadgs_LINK = $(CXXLINK)

quaderf_SOURCES = tquaderf.f
nodist_quaderf_SOURCES = second.f
quaderf_LINK = $(CXXLINK)

quadts_SOURCES = tquadts.f
nodist_quadts_SOURCES = second.f
quadts_LINK = $(CXXLINK)

roots_SOURCES = roots.f
nodist_roots_SOURCES = second.f
roots_LINK = $(CXXLINK)

demo: $(DEMO)

CLEANFILES = $(MODULES) $(DEMO)

endif
