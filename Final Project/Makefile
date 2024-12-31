# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -I src/utils -I src/pp_form -I src/b_form -I src/json

# Targets
all: test_task_a test_task_b test_task_c test_task_d test_task_e test_task_json test_task_curve_fitting test_task_linear_spline test_task_test

# Target for Task A
test_task_a: src/pp_form/linear_spline.cpp src/pp_form/cubic_spline.cpp src/utils/spline_utils.cpp src/taskA/test_task_a.cpp
	$(CXX) $(CXXFLAGS) -o test_task_a \
	src/pp_form/linear_spline.cpp \
	src/pp_form/cubic_spline.cpp \
	src/utils/spline_utils.cpp \
	src/taskA/test_task_a.cpp

# Target for Task B
test_task_b: src/pp_form/cubic_spline.cpp src/b_form/b_spline.cpp src/utils/spline_utils.cpp src/taskB/test_task_b.cpp
	$(CXX) $(CXXFLAGS) -o test_task_b \
	src/pp_form/cubic_spline.cpp \
	src/b_form/b_spline.cpp \
	src/utils/spline_utils.cpp \
	src/taskB/test_task_b.cpp

# Target for Task C
test_task_c: src/pp_form/cubic_spline.cpp src/b_form/b_spline.cpp src/utils/spline_utils.cpp src/taskC/test_task_c.cpp
	$(CXX) $(CXXFLAGS) -o test_task_c \
	src/pp_form/cubic_spline.cpp \
	src/b_form/b_spline.cpp \
	src/utils/spline_utils.cpp \
	src/taskC/test_task_c.cpp

# Target for Task D
test_task_d: src/pp_form/cubic_spline.cpp src/b_form/b_spline.cpp src/utils/spline_utils.cpp src/taskD/test_task_d.cpp
	$(CXX) $(CXXFLAGS) -o test_task_d \
	src/pp_form/cubic_spline.cpp \
	src/b_form/b_spline.cpp \
	src/utils/spline_utils.cpp \
	src/taskD/test_task_d.cpp

# Target for Task E
test_task_e: src/pp_form/cubic_spline.cpp src/b_form/b_spline.cpp src/utils/spline_utils.cpp src/taskE/test_task_e.cpp
	$(CXX) $(CXXFLAGS) -o test_task_e \
	src/pp_form/cubic_spline.cpp \
	src/b_form/b_spline.cpp \
	src/utils/spline_utils.cpp \
	src/taskE/test_task_e.cpp

# Target for JSON Task
test_task_json: src/json/json_config_handler.cpp src/taskJson/test_task_json.cpp
	$(CXX) $(CXXFLAGS) -o test_task_json \
	src/json/json_config_handler.cpp \
	src/taskJson/test_task_json.cpp

# Target for Linear Spline Task
test_task_check: src/pp_form/polynomial.cpp src/pp_form/piecewise_polynomial.cpp src/utils/math_function.cpp src/pp_form/pp_spline.cpp src/b_form/b_spline.cpp src/taskCheck/test_task_check.cpp
	$(CXX) $(CXXFLAGS) -o test_task_check \
	src/pp_form/polynomial.cpp \
	src/pp_form/piecewise_polynomial.cpp \
	src/utils/math_function.cpp \
	src/pp_form/pp_spline.cpp \
	src/b_form/b_spline.cpp \
	src/taskCheck/test_task_check.cpp

# Target for Curve Fitting Task
test_task_curve_fitting: src/pp_form/polynomial.cpp src/pp_form/piecewise_polynomial.cpp src/utils/math_function.cpp src/pp_form/pp_spline.cpp src/curve_fitting/curve_fitting.cpp src/taskCurveFitting/test_curve_fitting.cpp
	g++ -std=c++17 -I src/utils -I src/pp_form -I src/b_form -I src/json -o test_task_curve_fitting \
	src/pp_form/polynomial.cpp \
	src/pp_form/piecewise_polynomial.cpp \
	src/utils/math_function.cpp \
	src/pp_form/pp_spline.cpp \
	src/curve_fitting/curve_fitting.cpp \
	src/taskCurveFitting/test_curve_fitting.cpp

# Target for the new test
test_task_test: src/pp_form/polynomial.cpp src/pp_form/piecewise_polynomial.cpp src/utils/math_function.cpp src/pp_form/pp_spline.cpp src/b_form/b_spline.cpp src/taskTest/test_task_test.cpp
	$(CXX) $(CXXFLAGS) -o test_task_test \
	src/pp_form/polynomial.cpp \
	src/pp_form/piecewise_polynomial.cpp \
	src/utils/math_function.cpp \
	src/pp_form/pp_spline.cpp \
	src/b_form/b_spline.cpp \
	src/taskTest/test_task_test.cpp

# Clean target
clean:
	rm -f test_task_a test_task_b test_task_c test_task_d test_task_e test_task_json test_task_curve_fitting test_task_linear_spline test_task_test *.txt *.csv
