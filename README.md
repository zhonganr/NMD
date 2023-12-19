# Nonlinear Magnetoelastic Dynamics 
 
Advanced software for static and dynamic micromagnetic simulations, focusing on magnetoelastic interaction.

### Developers:
  Xiaoming LAN & Anruo ZHONG

### Main features:
  * Platform independent
  * CPU parallel computing
  * Minimal dependency on external libraries
  * Real-time control of the parameters
  * GUI based on Qt and OpenGL

### Dependencies:
  * gcc >=8.1.0 (possibly can be compiled also with earlier versions of gcc)
  * OpenGL >= 4.6
  * Qt >= 5.15

### Building from source:
  * Install a C compiler
    - Windows: [LLVM] (https://releases.llvm.org/)
  * Install [Qt](https://www.qt.io/download-open-source)
  * Compile with 
    - Windows: `g++ -std=c++17 -Wall -Wextra -g src/main.cpp src/moc_GUI.cpp src/GUI.cpp src/OPGL.cpp src/BoundaryShape.cpp src/Elastic_Field.cpp src/Energy.cpp src/Initial_State.cpp src/InOut.cpp src/Modify_MC.cpp src/Modify_NL.cpp src/Modify_Target.cpp src/GNEB.cpp src/LLG.cpp src/glad.c src/Array.cpp -I./include -I"${PATH_LLVM}\LLVM\include" -I"${PATH_LLVM}\LLVM\include\QtGui" -I"${PATH_LLVM}\LLVM\include\QtWidgets" -I"${PATH_LLVM}\LLVM\include\QtCore" -L./lib -L"${PATH_LLVM}\LLVM\lib" glfw3.dll -lfreetype -lz -lQt5Core -lQt5Gui -lQt5Widgets -o bin/NMD.exe` or with Makefile: `mingw32-make run`
