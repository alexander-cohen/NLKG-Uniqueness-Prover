CONFIG_FILE = LinuxWithProfil
INSTALL_DIR = $(HOME)/VNODELP_INSTALL

include $(INSTALL_DIR)/vnodelp/config/$(CONFIG_FILE)

CXXFLAGS += -I$(INSTALL_DIR)/vnodelp/include # compiler flags
CXXFLAGS += -I$(INSTALL_DIR)/vnodelp/FADBAD++  
CXXFLAGS += -std=c++11 -g
LDFLAGS  += -L$(INSTALL_DIR)/vnodelp/lib    # library flags
LIBS     = -lvnode $(I_LIBS) $(LAPACK_LIB) \
    $(BLAS_LIB) $(GPP_LIBS) 

OBJS = nlkg_uniqueness_prover.o

all: solver

solver: $(OBJS) nlkg_helper.h
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	@-$(RM) *.o *.out core.* solver
