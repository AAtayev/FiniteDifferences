# input files
# TEX=report.tex
# PLOTSCRIPTS=plotscript1.gpl plotscript2.gpl

# generated files
# REPORT=report.pdf
RESULTS=Test.dat
PROGRAM=FiniteDifference
OBJS=main.o FiniteDifference.o SparseMatrix.o

# additional variables
CPPFLAGS=-std=c++11

all: $(RESULTS)

$(RESULTS): $(PROGRAM)
	./$(PROGRAM)

$(PROGRAM): $(OBJS)
	g++ $(CPPFLAGS) $(OBJS) -o $(PROGRAM)

$(OBJS): %.o: %.cc
	g++ -Ofast -Wall -Wfatal-errors $(CPPFLAGS) -c $^ -o  $@

clean:
	rm -rf $(OBJS) $(PROGRAM) $(RESULTS)
