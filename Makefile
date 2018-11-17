# input files
# TEX=report.tex
PLOTSCRIPTS=plotscript1.gpl plotscript2.gpl plotscript3.gpl plotscript4.gpl plotscript5.gpl

# generated files
# REPORT=report.pdf
RESULTS1=Error_Pe_0.000000.dat Error_Pe_0.000500.dat Error_Pe_1.000000.dat Error_Pe_10.500000.dat Solution_Pe_0.000000_J_9.dat Solution_Pe_0.000000_J_19.dat Solution_Pe_0.000000_J_39.dat Solution_Pe_0.000000_J_79.dat Solution_Pe_0.000000_J_159.dat Solution_Pe_0.000000_J_319.dat Solution_Pe_0.000000_J_639.dat Solution_Pe_0.000500_J_9.dat Solution_Pe_0.000500_J_19.dat Solution_Pe_0.000500_J_39.dat Solution_Pe_0.000500_J_79.dat Solution_Pe_0.000500_J_159.dat Solution_Pe_0.000500_J_319.dat Solution_Pe_0.000500_J_639.dat Solution_Pe_1.000000_J_9.dat Solution_Pe_1.000000_J_19.dat Solution_Pe_1.000000_J_39.dat Solution_Pe_1.000000_J_79.dat Solution_Pe_1.000000_J_159.dat Solution_Pe_1.000000_J_319.dat Solution_Pe_1.000000_J_639.dat Solution_Pe_10.500000_J_9.dat Solution_Pe_10.500000_J_19.dat Solution_Pe_10.500000_J_39.dat Solution_Pe_10.500000_J_79.dat Solution_Pe_10.500000_J_159.dat Solution_Pe_10.500000_J_319.dat Solution_Pe_10.500000_J_639.dat
RESULTS=$(wildcard *.dat)
PROGRAM=FiniteDifference
OBJS=main.o FiniteDifference.o SparseMatrix.o
PLOTS=Error_plots.pdf Solution_plots_Pe_0.pdf Solution_plots_Pe_0.0005.pdf Solution_plots_Pe_1.pdf Solution_plots_Pe_10.5.pdf

# additional variables
CPPFLAGS=-std=c++11

all: $(PLOTS)

$(PLOTS): $(RESULTS1) $(RESULTS) $(PLOTSCRIPTS)
	gnuplot plotscript1.gpl
	gnuplot plotscript2.gpl
	gnuplot plotscript3.gpl
	gnuplot plotscript4.gpl
	gnuplot plotscript5.gpl

$(RESULTS1): $(PROGRAM)
	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -Ofast -Wall -Wfatal-errors $(CPPFLAGS) -c $^ -o  $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(RESULTS)
