CXX_SOURCES = A.cpp B.cpp C.cpp D.cpp E.cpp F.cpp
CXX_EXECUTABLES = $(CXX_SOURCES:.cpp=)

PDF_TARGET = HW2_Programing.pdf
PYTHON_SCRIPT = zuotu.py

all: run_cpp run_python pdf

$(CXX_EXECUTABLES): $(CXX_SOURCES)
	g++ -o $@ $(@:=.cpp)
run_cpp: $(CXX_EXECUTABLES)
	@for exe in $(CXX_EXECUTABLES); do \
		./$$exe; \
	done

run_python:
	@echo "Running $(PYTHON_SCRIPT)..."
	python3 $(PYTHON_SCRIPT)

pdf: $(PDF_TARGET)
$(PDF_TARGET): HW2_Programing.tex
	xelatex HW2_Programing.tex
clean:
	rm -f $(PDF_TARGET) $(CXX_EXECUTABLES) *.aux *.log
.PHONY: all clean run_cpp run_python pdf
