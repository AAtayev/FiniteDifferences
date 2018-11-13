# input files
# TEX=report.tex
PLOTSCRIPTS=plotscript1.gpl

# generated files
# REPORT=report.pdf
RESULTS1=Error_Pe_0.000500.dat Error_Pe_1.000000.dat Error_Pe_10.500000.dat
RESULTS=$(wildcard *.dat)
PROGRAM=FiniteDifference
OBJS=main.o FiniteDifference.o SparseMatrix.o
PLOTS=Plots.pdf

# additional variables
CPPFLAGS=-std=c++11

all: $(PLOTS)

$(PLOTS): $(RESULTS1) $(PLOTSCRIPTS)
	gnuplot plotscript1.gpl

$(RESULTS1): $(PROGRAM)
	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -Ofast -Wall -Wfatal-errors $(CPPFLAGS) -c $^ -o  $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(RESULTS)
