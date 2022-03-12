
cmake_minimum_required(VERSION 3.10.0)
cmake_policy(SET CMP0020 NEW)
cmake_policy(SET CMP0076 NEW)
set(CMAKE_SOURCE_DIR src)

set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION ".")
include(InstallRequiredSystemLibraries)

set(PROJECT_NAME "ConeSegmentation")
set(PROJECT_VER "1.1.0")
set(PROJECT_BUILD "03")
PROJECT(${PROJECT_NAME})

string(TIMESTAMP PROJECT_DATE "%Y-%m-%d")
message("Configuring ${PROJECT_NAME} ver. ${PROJECT_VER}${PROJECT_BUILD} (${PROJECT_DATE})")

file(WRITE src/version.h
     "/* Auto-generated file. Don't modify! */\n"
     "const char *__version__ = \"${PROJECT_VER}${PROJECT_BUILD} (${PROJECT_DATE})\";\n"
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Give preference to static libraries when possible
set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

### Qt
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets PrintSupport)

### Vtk
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

### Itk
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

### Main project
file(GLOB UI_FILES src/*.ui)
file(GLOB QT_WRAP src/*.h)
file(GLOB CXX_FILES src/*.cxx)
file(GLOB CPP_FILES src/*.cpp)
file(GLOB TXX_FILES src/*.tpp)

message("QT5: ${Qt5_DIR}")
message("VTK: ${VTK_DIR}")
message("ITK: ${ITK_DIR}")

set(VTK_VER "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
message("VTK_VER: ${VTK_VER}")

add_executable(${PROJECT_NAME} WIN32 MACOSX_BUNDLE ${CXX_FILES} ${TXX_FILES} ${CPP_FILES} Icons/qresources.qrc)
 
file(COPY "Help" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

target_link_libraries(${PROJECT_NAME}
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  Qt5::Core
  Qt5::Gui
  Qt5::Widgets
)

if(WIN32)
target_sources(${PROJECT_NAME} PUBLIC src/ConeSegmentation.rc)

target_link_libraries(${PROJECT_NAME}
  Qt5::WinMain
)

set(VTK_DLL ${VTK_DIR}/bin)
set(QT5_DLL ${Qt5_DIR}/../../../bin)

macro(cp_dlls0 APP DLL)
    # find the release *.dll file
    get_target_property(Qt5_${DLL}Location Qt5::${DLL} LOCATION)
    # find the debug *d.dll file
    get_target_property(Qt5_${DLL}LocationDebug Qt5::${DLL} IMPORTED_LOCATION_DEBUG)
    add_custom_command(TARGET ${APP} POST_BUILD
       COMMAND ${CMAKE_COMMAND} -E copy_if_different $<$<CONFIG:Debug>:${Qt5_${DLL}LocationDebug}> $<$<NOT:$<CONFIG:Debug>>:${Qt5_${DLL}Location}> $<TARGET_FILE_DIR:${APP}>)
        #COMMAND ${CMAKE_COMMAND} -E copy_if_different $<IF $<CONFIG:Debug>, ${Qt5_${DLL}LocationDebug}, ${Qt5_${DLL}Location}> $<TARGET_FILE_DIR:${APP}>)
endmacro()

macro(cp_vtk APP DLL)
    add_custom_command(TARGET ${APP} POST_BUILD
       COMMAND ${CMAKE_COMMAND} -E copy_if_different ${VTK_DLL}/$<CONFIG>/${DLL} $<TARGET_FILE_DIR:${APP}>)
endmacro()

macro(cp_qt5 APP DLL)
    add_custom_command(TARGET ${APP} POST_BUILD
       COMMAND ${CMAKE_COMMAND} -E copy_if_different $<$<CONFIG:Debug>:${DLL}d.dll> $<$<NOT:$<CONFIG:Debug>>:${DLL}.dll> $<TARGET_FILE_DIR:${APP}>)
endmacro()

macro(cp_qt5_other APP DLL)
    add_custom_command(TARGET ${APP} POST_BUILD
       COMMAND ${CMAKE_COMMAND} -E copy_if_different ${DLL} $<TARGET_FILE_DIR:${APP}>)
endmacro()

# copy all xxd.dll from QT5 depending on msvc debug or release build
list(APPEND QT5_DLLS
    "${QT5_DLL}/Qt5Core"
    "${QT5_DLL}/Qt5Gui"
    "${QT5_DLL}/Qt5Widgets"
    )
foreach( dll ${QT5_DLLS} )
    cp_qt5(${PROJECT_NAME} ${dll})
    install(FILES ${dll}.dll DESTINATION . CONFIGURATIONS Release)
endforeach(dll)

# copy all vtk dll depending on msvc debug or release build
#
file(GLOB VTK_DLLS "${VTK_DLL}/Release/*-${VTK_VER}.dll")

foreach( dll_path ${VTK_DLLS} )
	get_filename_component(dll ${dll_path} NAME)
    cp_vtk(${PROJECT_NAME} ${dll})
    install(FILES ${VTK_DLL}/Release/${dll} DESTINATION . CONFIGURATIONS Release)
endforeach(dll_path)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${Qt5_DIR}/../../../plugins/platforms/qwindows.dll" "$<TARGET_FILE_DIR:${PROJECT_NAME}>/platforms/qwindows.dll")
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${Qt5_DIR}/../../../plugins/imageformats/qtiff.dll" "$<TARGET_FILE_DIR:${PROJECT_NAME}>/imageformats/qtiff.dll")
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${Qt5_DIR}/../../../translations/qt_en.qm" "$<TARGET_FILE_DIR:${PROJECT_NAME}>/translations/qt_en.qm")

install(FILES "${Qt5_DIR}/../../../plugins/platforms/qwindows.dll" DESTINATION "./platforms" CONFIGURATIONS Release)
#install(FILES "${Qt5_DIR}/../../../plugins/iconengines/qsvgicon.dll" DESTINATION "./iconengines")
install(FILES "${Qt5_DIR}/../../../plugins/imageformats/qtiff.dll" DESTINATION "./imageformats" CONFIGURATIONS Release)
install(FILES "${Qt5_DIR}/../../../translations/qt_en.qm" DESTINATION "./translations" CONFIGURATIONS Release)
install(DIRECTORY "./Help" DESTINATION ".")
install(TARGETS ${PROJECT_NAME} DESTINATION . CONFIGURATIONS Release)

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "NSIS") # this needs NSIS installed, and available
#set(CPACK_SOURCE_GENERATOR "ZIP")
set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
#set(CPACK_BUILD_SOURCE_DIRS "src;${CMAKE_BINARY_DIR}/${Debug_Release}")
set(CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE 1)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_NAME}")
set(CPACK_NSIS_MUI_UNIICON "${CMAKE_CURRENT_SOURCE_DIR}/Icons/SegmentIcon256x256.ico")
set (CPACK_RESOURCE_FILE_LICENSE
     "${CMAKE_CURRENT_SOURCE_DIR}/Help/License.txt")
set (CPACK_PACKAGE_VENDOR "National Eye Institute")
set (CPACK_PACKAGE_VERSION ${PROJECT_VER})
set (CPACK_PACKAGE_EXECUTABLES "${PROJECT_NAME}" "Cone Segmentation")
set (CPACK_NSIS_EXECUTABLES_DIRECTORY ".")
set (CPACK_CREATE_DESKTOP_LINKS "${PROJECT_NAME}")
set (CPACK_NSIS_MODIFY_PATH "ON")
set(CPACK_SOURCE_IGNORE_FILES "Debug")
include (CPack)

endif()
if (APPLE)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9" CACHE STRING "Minimum MacOS X version required for the project" FORCE)

# set(CPACK_DMG_DISABLE_APPLICATIONS_SYMLINK ON)
set(MACOSX_BUNDLE_ICON_FILE SegmentIcon256x256.icns)
set(APP_ICON "${CMAKE_CURRENT_SOURCE_DIR}/Icons/SegmentIcon256x256.icns")
set_source_files_properties(${APP_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
target_sources(${PROJECT_NAME} PUBLIC ${APP_ICON})

set(prefix "${PROJECT_NAME}.app/Contents")
set(INSTALL_RUNTIME_DIR "${prefix}/MacOS")
set(INSTALL_CMAKE_DIR "${prefix}/Resources")

# based on code from CMake's QtDialog/CMakeLists.txt
macro(install_qt5_plugin _qt_plugin_name _qt_plugins_var _prefix)
    get_target_property(_qt_plugin_path "${_qt_plugin_name}" LOCATION)
    if(EXISTS "${_qt_plugin_path}")
        get_filename_component(_qt_plugin_file "${_qt_plugin_path}" NAME)
        get_filename_component(_qt_plugin_type "${_qt_plugin_path}" PATH)
        get_filename_component(_qt_plugin_type "${_qt_plugin_type}" NAME)
        set(_qt_plugin_dest "${_prefix}/PlugIns/${_qt_plugin_type}")
        install(FILES "${_qt_plugin_path}"
            DESTINATION "${_qt_plugin_dest}")
        set(${_qt_plugins_var}
            "${${_qt_plugins_var}};\$ENV{DEST_DIR}\${CMAKE_INSTALL_PREFIX}/${_qt_plugin_dest}/${_qt_plugin_file}")
    else()
        message(FATAL_ERROR "QT plugin ${_qt_plugin_name} not found")
    endif()
endmacro()

install_qt5_plugin("Qt5::QCocoaIntegrationPlugin" QT_PLUGINS ${prefix})
install_qt5_plugin("Qt5::QTiffPlugin" QT_PLUGINS ${prefix})

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/Help" DESTINATION ${prefix})

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
    "[Paths]\nPlugins = ${_qt_plugin_dir}\n")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
    DESTINATION "${INSTALL_CMAKE_DIR}")

# Destination paths below are relative to ${CMAKE_INSTALL_PREFIX}
install(TARGETS ${PROJECT_NAME}
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION ${INSTALL_RUNTIME_DIR} COMPONENT Runtime
    )

# Note Mac specific extension .app
set(APPS "\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}.app")

# Directories to look for dependencies
set(DIRS "${CMAKE_BINARY_DIR}")

# Path used for searching by FIND_XXX(), with appropriate suffixes added
if(CMAKE_PREFIX_PATH)
    foreach(dir ${CMAKE_PREFIX_PATH})
        list(APPEND DIRS "${dir}/bin" "${dir}/lib")
    endforeach()
endif()

# Append Qt's lib folder which is two levels above Qt5Widgets_DIR
list(APPEND DIRS "${Qt5Widgets_DIR}/../..")

include(InstallRequiredSystemLibraries)

message(STATUS "APPS: ${APPS}")
message(STATUS "QT_PLUGINS: ${QT_PLUGINS}")
message(STATUS "DIRS: ${DIRS}")

install(CODE "include(BundleUtilities)
    fixup_bundle(\"${APPS}\" \"${QT_PLUGINS}\" \"${DIRS}\")")

set(CPACK_PACKAGE_VERSION ${PROJECT_VER})
# set(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/Icons/SegmentIcon256x256.icns")
set (CPACK_RESOURCE_FILE_LICENSE
     "${CMAKE_CURRENT_SOURCE_DIR}/Help/License.txt")
set (CPACK_PACKAGE_VENDOR "National Eye Institute")

set(CPACK_DMG_BACKGROUND_IMAGE "${CMAKE_CURRENT_SOURCE_DIR}/Packaging/DMGBackground.tif")
set(CPACK_DMG_DS_STORE_SETUP_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/Packaging/DMGSetup.scpt")

set(CPACK_GENERATOR "DragNDrop")
include(CPack)
endif()
