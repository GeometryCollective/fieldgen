#######################################################################################################

# Specify library locations here (add or remove "#" marks to comment/uncomment lines for your platform)

# Mac OS X
DDG_INCLUDE_PATH      = -I/usr/local/include
DDG_LIBRARY_PATH      = -L/usr/local/lib
DDG_BLAS_LIBS         = -framework Accelerate
DDG_SUITESPARSE_LIBS  = -lspqr -lumfpack -lcholmod -lmetis -lcolamd -lccolamd -lcamd -lamd -lm -lsuitesparseconfig
DDG_OPENGL_LIBS       = -framework OpenGL -framework GLUT

# # Linux
# DDG_INCLUDE_PATH      =
# DDG_LIBRARY_PATH      =
# DDG_BLAS_LIBS         = -llapack -lblas -lgfortran
# DDG_SUITESPARSE_LIBS  = -lspqr -lcholmod -lmetis -lcolamd -lccolamd -lcamd -lamd -lm
# DDG_OPENGL_LIBS       = -lglut -lGL -lGLU -lX11

# # Windows / Cygwin
# DDG_INCLUDE_PATH      = -I/usr/include/opengl -I/usr/include/suitesparse
# DDG_LIBRARY_PATH      = -L/usr/lib/w32api -L/usr/lib/suitesparse
# DDG_BLAS_LIBS         = -llapack -lblas
# DDG_SUITESPARSE_LIBS  = -lspqr -lcholmod -lcolamd -lccolamd -lcamd -lamd -lm
# DDG_OPENGL_LIBS       = -lglut32 -lglu32 -lopengl32

#######################################################################################################

TARGET = fieldviz
CLTARGET = fieldgen
CC = g++
LD = g++
CFLAGS = -g -std=c++11 -Wall -Wno-deprecated -Werror -Wno-error=deprecated-declarations -ansi -pedantic  $(DDG_INCLUDE_PATH) -I./include -I./src
LFLAGS = -g -Wall -Wno-deprecated -Werror -pedantic $(DDG_LIBRARY_PATH)
LIBS = $(DDG_OPENGL_LIBS) $(DDG_SUITESPARSE_LIBS) $(DDG_BLAS_LIBS)
CLLIBS = $(DDG_SUITESPARSE_LIBS) $(DDG_BLAS_LIBS)

## !! Do not edit below this line -- dependencies can be updated by running ./update ##################

OBJS = obj/Camera.o obj/Complex.o obj/DenseMatrix.o obj/Edge.o obj/Face.o obj/HalfEdge.o obj/Image.o obj/KVecDir.o obj/LinearContext.o obj/Mesh.o obj/MeshIO.o obj/Quaternion.o obj/Real.o obj/SectionIntegrals.o obj/Shader.o obj/SparseMatrix.o obj/Vector.o obj/Vertex.o obj/Viewer.o obj/main.o
CLOBJS = obj/Complex.o obj/DenseMatrix.o obj/Edge.o obj/Face.o obj/HalfEdge.o obj/KVecDir.o obj/LinearContext.o obj/Mesh.o obj/MeshIO.o obj/Quaternion.o obj/Real.o obj/SectionIntegrals.o obj/SparseMatrix.o obj/Vector.o obj/Vertex.o obj/commandline.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) $(LIBS) -o $(TARGET)

commandline: $(CLOBJS)
	$(LD) $(LFLAGS) $(CLOBJS) $(CLLIBS) -o $(CLTARGET)

obj/Camera.o: src/Camera.cpp include/Camera.h include/Quaternion.h include/Vector.h 
	$(CC) $(CFLAGS) -c src/Camera.cpp -o obj/Camera.o

obj/Complex.o: src/Complex.cpp include/Complex.h 
	$(CC) $(CFLAGS) -c src/Complex.cpp -o obj/Complex.o

obj/DenseMatrix.o: src/DenseMatrix.cpp include/DenseMatrix.h include/Types.h src/DenseMatrix.inl include/DenseMatrix.h include/LinearContext.h include/Quaternion.h include/Vector.h include/SparseMatrix.h include/Complex.h src/SparseMatrix.inl include/Real.h include/Complex.h include/Utility.h 
	$(CC) $(CFLAGS) -c src/DenseMatrix.cpp -o obj/DenseMatrix.o

obj/Edge.o: src/Edge.cpp include/Edge.h include/Types.h include/Quaternion.h include/Vector.h include/Complex.h include/Mesh.h include/HalfEdge.h include/Vertex.h include/Edge.h include/Face.h include/AliasTable.h
	$(CC) $(CFLAGS) -c src/Edge.cpp -o obj/Edge.o

obj/Face.o: src/Face.cpp include/Face.h include/Types.h include/Complex.h include/Vector.h include/Mesh.h include/HalfEdge.h include/Quaternion.h include/Vertex.h include/Edge.h include/Face.h include/AliasTable.h include/Vector.h include/Utility.h 
	$(CC) $(CFLAGS) -c src/Face.cpp -o obj/Face.o

obj/HalfEdge.o: src/HalfEdge.cpp include/HalfEdge.h include/Vector.h include/Types.h include/Complex.h include/Quaternion.h include/Edge.h include/Mesh.h include/HalfEdge.h include/Vertex.h include/Edge.h include/Face.h include/AliasTable.h include/Vertex.h include/Vector.h 
	$(CC) $(CFLAGS) -c src/HalfEdge.cpp -o obj/HalfEdge.o

obj/Image.o: src/Image.cpp include/Image.h 
	$(CC) $(CFLAGS) -c src/Image.cpp -o obj/Image.o

obj/KVecDir.o: src/KVecDir.cpp include/Utility.h include/Complex.h include/Complex.h include/Mesh.h include/HalfEdge.h include/Vector.h include/Types.h include/Quaternion.h include/Vertex.h include/Edge.h include/Face.h include/AliasTable.h include/SparseMatrix.h src/SparseMatrix.inl include/Real.h include/Complex.h include/SparseMatrix.h include/DenseMatrix.h src/DenseMatrix.inl include/LinearContext.h include/Quaternion.h include/Utility.h include/SectionIntegrals.h 
	$(CC) $(CFLAGS) -c src/KVecDir.cpp -o obj/KVecDir.o

obj/LinearContext.o: src/LinearContext.cpp include/LinearContext.h 
	$(CC) $(CFLAGS) -c src/LinearContext.cpp -o obj/LinearContext.o

obj/Mesh.o: src/Mesh.cpp include/Mesh.h include/HalfEdge.h include/Vector.h include/Types.h include/Complex.h include/Quaternion.h include/Vertex.h include/Edge.h include/Face.h include/AliasTable.h include/MeshIO.h 
	$(CC) $(CFLAGS) -c src/Mesh.cpp -o obj/Mesh.o

obj/MeshIO.o: src/MeshIO.cpp include/MeshIO.h include/Vector.h include/Mesh.h include/HalfEdge.h include/Types.h include/Complex.h include/Quaternion.h include/Vertex.h include/Edge.h include/Face.h include/AliasTable.h 
	$(CC) $(CFLAGS) -c src/MeshIO.cpp -o obj/MeshIO.o

obj/Quaternion.o: src/Quaternion.cpp include/Quaternion.h include/Vector.h include/Complex.h 
	$(CC) $(CFLAGS) -c src/Quaternion.cpp -o obj/Quaternion.o

obj/Real.o: src/Real.cpp include/Real.h 
	$(CC) $(CFLAGS) -c src/Real.cpp -o obj/Real.o

obj/SectionIntegrals.o: src/SectionIntegrals.cpp include/Complex.h include/SectionIntegrals.h include/Complex.h 
	$(CC) $(CFLAGS) -c src/SectionIntegrals.cpp -o obj/SectionIntegrals.o

obj/Shader.o: src/Shader.cpp include/Shader.h 
	$(CC) $(CFLAGS) -c src/Shader.cpp -o obj/Shader.o

obj/SparseMatrix.o: src/SparseMatrix.cpp include/SparseMatrix.h include/Types.h include/Complex.h src/SparseMatrix.inl include/Real.h include/Complex.h include/SparseMatrix.h include/DenseMatrix.h src/DenseMatrix.inl include/LinearContext.h include/Quaternion.h include/Vector.h include/Utility.h 
	$(CC) $(CFLAGS) -c src/SparseMatrix.cpp -o obj/SparseMatrix.o

obj/Vector.o: src/Vector.cpp include/Vector.h include/Complex.h 
	$(CC) $(CFLAGS) -c src/Vector.cpp -o obj/Vector.o

obj/Vertex.o: src/Vertex.cpp include/Vertex.h include/Vector.h include/Types.h include/Complex.h include/HalfEdge.h include/Quaternion.h include/Mesh.h include/Vertex.h include/Edge.h include/Face.h include/AliasTable.h include/HalfEdge.h include/Quaternion.h 
	$(CC) $(CFLAGS) -c src/Vertex.cpp -o obj/Vertex.o

obj/Viewer.o: src/Viewer.cpp include/Viewer.h include/Mesh.h include/HalfEdge.h include/Vector.h include/Types.h include/Complex.h include/Quaternion.h include/Vertex.h include/Edge.h include/Face.h include/AliasTable.h include/Camera.h include/Shader.h include/Image.h include/Utility.h include/SectionIntegrals.h 
	$(CC) $(CFLAGS) -c src/Viewer.cpp -o obj/Viewer.o

obj/main.o: src/main.cpp include/Viewer.h include/Mesh.h include/HalfEdge.h include/Vector.h include/Types.h include/Complex.h include/Quaternion.h include/Vertex.h include/Edge.h include/Face.h include/AliasTable.h include/Camera.h include/Shader.h include/DenseMatrix.h src/DenseMatrix.inl include/DenseMatrix.h include/LinearContext.h include/Quaternion.h include/SparseMatrix.h src/SparseMatrix.inl include/Real.h include/Complex.h include/Utility.h 
	$(CC) $(CFLAGS) -c src/main.cpp -o obj/main.o

obj/commandline.o: src/main.cpp include/Mesh.h include/HalfEdge.h include/Vector.h include/Types.h include/Complex.h include/Quaternion.h include/Vertex.h include/Edge.h include/Face.h include/Camera.h include/DenseMatrix.h src/DenseMatrix.inl include/DenseMatrix.h include/LinearContext.h include/Quaternion.h include/SparseMatrix.h src/SparseMatrix.inl include/Real.h include/Complex.h include/Utility.h 
	$(CC) $(CFLAGS) -c src/commandline.cpp -o obj/commandline.o


clean:
	rm -f $(OBJS)
	rm -f $(CLOBJS)
	rm -f $(TARGET)
	rm -f $(TARGET).exe
	rm -f $(CLTARGET)
	rm -f $(CLTARGET).exe

