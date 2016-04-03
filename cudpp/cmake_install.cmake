# Install script for directory: /Users/zen/Dropbox/DCS/Experiment/cudpp

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE RENAME "cudpp-config.cmake" FILES "/Users/zen/Dropbox/DCS/Experiment/eulercuda_cython/cudpp/lib/cudpp-config.cmake.install")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "/Users/zen/Dropbox/DCS/Experiment/eulercuda_cython/cudpp/lib/cudpp-config-version.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/zen/Dropbox/DCS/Experiment/eulercuda_cython/cudpp/src/cudpp/cmake_install.cmake")
  include("/Users/zen/Dropbox/DCS/Experiment/eulercuda_cython/cudpp/src/cudpp_hash/cmake_install.cmake")
  include("/Users/zen/Dropbox/DCS/Experiment/eulercuda_cython/cudpp/apps/cudpp_testrig/cmake_install.cmake")
  include("/Users/zen/Dropbox/DCS/Experiment/eulercuda_cython/cudpp/apps/cudpp_hash_testrig/cmake_install.cmake")
  include("/Users/zen/Dropbox/DCS/Experiment/eulercuda_cython/cudpp/apps/simpleCUDPP/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/zen/Dropbox/DCS/Experiment/eulercuda_cython/cudpp/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
