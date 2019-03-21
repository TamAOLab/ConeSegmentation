
Building for Windows:
   *NOTE* VTK is built and linked dynamically

1. Install MS Visual Studio 2017
2. Install CMake and NSIS
3. Download and install QT 5 (Checking "msvc2017 64-bit" under the latest 5.x version is sufficient)
4. Download and build VTK 8 (both Debug and Release) with the following CMake options checked:
   `VTK_Group_Qt`
   `VTK_LEGACY_SILENT` (Advanced)
5. Download and build ITK 4 (both Debug and Release) with the following CMake options checked:
   `Module_ITKV3Compatibility`
   `ITKV3_COMPATIBILITY`
   `Module_ITKVtkGlue` (Advanced)
6. Build cone-segmentation by running CMake first and providing the following options:
   `Qt5_DIR`:PATH - Qt cmake dir (containing Qt5Config.cmake), e.g. "C:/Qt/5.11.1/msvc2017_64/lib/cmake/Qt5/"
   `VTK_DIR`:PATH - VTK build dir (containing CMakeCache.txt), e.g. "C:/VTK/VTKBuild/"
   `ITK_DIR`:PATH - ITK build dir (containing CMakeCache.txt), e.g. "C:/ITK/ITKBuild/"
7. To make the distribution package with VisualStudio, switch to "Release" mode,
then build target PACKAGE. The result is "ConeSegmentation-<X.X.X>-win64.exe" in the build directory.

Building for MacOS:
   *NOTE* VTK is built and linked statically

1. Install Xcode
2. Install CMake
3. Download and install Qt 5 (Checking "macOS" under the latest 5.x version is sufficient)
4. Download and build VTK 8 (both Debug and Release) with the following CMake options checked:
   `VTK_Group_Qt`
   `VTK_LEGACY_SILENT` (Advanced)
and the `BUILD_SHARED_LIBS` option *unchecked* (unlike in the Windows build!)
5. Download and build ITK 4 (both Debug and Release) with the following CMake options checked:
   `Module_ITKV3Compatibility`
   `ITKV3_COMPATIBILITY`
   `Module_ITKVtkGlue` (Advanced)
6. Build cone-segmentation project by running CMake, selecting Xcode as compiler,
and providing the following options:
   `Qt5_DIR`:PATH - Qt cmake dir (containing Qt5Config.cmake), e.g. "/Users/Shared/Qt/5.12.0/clang_64/lib/cmake/Qt5/"
   `VTK_DIR`:PATH - VTK build dir (containing CMakeCache.txt), e.g. "/Users/Shared/VTK/build/"
   `ITK_DIR`:PATH - ITK build dir (containing CMakeCache.txt), e.g. "/Users/Shared/ITK/build/"
7. To make the distribution package with Xcode, switch ALL_BUILD, install and PACKAGE targets to Release mode,
then build target PACKAGE. The result is "ConeSegmentation-<X.X.X>-Darwin.dmg" in the build directory.
If Xcode asks for permission to access "Finder", answer "Yes". 

