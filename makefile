#---------------------------------------------------------------
# Directories Of the project
#---------------------------------------------------------------
SRCDIR=src
OBJDIR=obj
EXECDIR=exec
DATADIR=data
MODSDIR=mods
IMGDIR=images

EXE=$(EXECDIR)/HLLC

VPATH=$(SRC):$(OBJ):$(EXEC)

FFLAGS= -U -03
FC=gfortran $(FFLAGS)


OBJ=\
parameters.o\
globals.o\
hydrocore.o\
utilities.o\
init.o\
boundaries.o\
hllc.o\
solver.o

MAIN=$(OBJDIR)/main.o

OBJS=$(patsubst %.o, $(OBJDIR)/%.o, $(OBJ))
SRCS=$(patsubst $(OBJDIR)/%.o, $(SRCDIR)/%.f03, $(OBJS))
MODS=$(patsubst $(OBJDIR)/%.o, $(MODSDIR)/%.mod, $(OBJS))

.PHONY= all run rerun clean plot gif compilingMessage

all: compile run plot gif

rerun: all | clean

compile: compilingMessage $(EXE)

$(EXE): $(OBJS) $(MAIN) | $(MODS) $(EXECDIR) $(MODSDIR)
	@$(FC) -I$(MODSDIR) $^ -o $@

$(MAIN): $(patsubst $(OBJDIR)/%.o, $(SRCDIR)/%.f03, $(MAIN))
	@echo "Compiling " $^ " ..."
	@$(FC) -I$(MODSDIR) -c $< -o $@

$(MODSDIR)/%.mod $(OBJDIR)/%.o: $(SRCDIR)/%.f03 | $(OBJDIR) $(MODSDIR)
	@echo "Compiling " $^ " ..."
	@$(FC) -J$(MODSDIR) -c $< -o $@

$(OBJDIR) $(MODSDIR) $(EXECDIR) $(DATADIR) $(IMGDIR):
	@mkdir $@

compilingMessage:
	@echo "=========================================================="
	@echo "Compiling the source files"
	@echo "=========================================================="

run: $(EXE) | $(DATADIR)
	@echo "=========================================================="
	@echo "Excecuting the program"
	@echo "=========================================================="
ifneq ("$(wildcard $(DATADIR)/*)","")
	rm $(DATADIR)/*
endif
ifneq ("$(wildcard $(IMGDIR)/*)","")
	rm $(IMGDIR)/*
endif
	@$<
	@mv *.dat $(DATADIR)

plot: $(DATADIR) $(IMGDIR)
	@echo "=========================================================="
	@echo "Generating images of the data"
	@echo "=========================================================="
	@for file in $(DATADIR)/*.dat; do \
			python $(SRCDIR)/plotter.py $$file ;\
			echo "ploted: " $$file ;\
			mv $(DATADIR)/*.png $(IMGDIR);\
	 done

gif:
	@echo "=========================================================="
	@echo "Generating animation of the data"
	@echo "=========================================================="
	@convert -delay 15 -loop 0 $(IMGDIR)/*.png simulation.gif
	@echo "Done."

clean: $(OBJDIR) $(MODSDIR) $(EXECDIR) $(DATADIR) $(IMGDIR)
	@echo "=========================================================="
	@echo "Cleaning the project"
	@echo "=========================================================="
	@rm -R $(OBJDIR) $(MODSDIR) $(EXECDIR) $(DATADIR) $(IMGDIR)
