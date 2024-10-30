CXX_SOURCES = A.cpp B.cpp C.cpp D.cpp E.cpp F.cpp
CXX_EXECUTABLES = $(CXX_SOURCES:.cpp=)

all: $(PDF_TARGET) run_cpp
xelatex HW2_Programing.tex


$(CXX_EXECUTABLES): $(CXX_SOURCES)
	g++ -o $@ $(@:=.cpp)

run_cpp: $(CXX_EXECUTABLES)
	@for exe in $(CXX_EXECUTABLES); do \
		./$$exe; \
	done


clean:
	rm -f $(PDF_TARGET) $(CXX_EXECUTABLES) *.aux *.log
.PHONY: all clean run_cpp
