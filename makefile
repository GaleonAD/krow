# location of the Python header files
 
#PYTHON_VERSION = 2.7
#PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
 
# location of the Boost Python include files and library
 
#BOOST_INC = /usr/include
#BOOST_LIB = /usr/lib

CODE = main.cc spherical_grid.cc 
TARGET = play
 
$(TARGET): $(CODE)
	g++ -std=c++11 $(CODE) -o $(TARGET)
