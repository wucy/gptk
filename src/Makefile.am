GPTK_RELEASE=1:0:0
GPTK_LIBRARY_VERSION=0:2:0

# The final library
lib_LTLIBRARIES = libgptk.la

# The sources (none, as they are all in subdirectories and will be
# compiled in the sub-libraries)
libgptk_la_SOURCES = 

# The sub-libraries
libgptk_la_LIBADD = gaussian_processes/libgp.la \
                    covariance_functions/libcovf.la \
                    optimisation/liboptim.la \
                    likelihood_models/liblikmod.la \
                    parameter_transforms/libtransforms.la \
                    itppext/libitppext.la \
                    io/libio.la \
                    plotting/libplot.la

# Some package information (version)
libgptk_la_LDFLAGS= -version-info $(GPTK_LIBRARY_VERSION) -release $(GPTK_RELEASE)

# Also run make in the following sub-directories 
SUBDIRS = gaussian_processes \
          covariance_functions \
          parameter_transforms \
          likelihood_models \
          optimisation \
          io \
          itppext \
          plotting
           
