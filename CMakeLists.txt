cmake_minimum_required(VERSION 3.10)

### --- --- --- --- --- --- Version --- --- --- --- --- --- ###
set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
set(VERSION_NO ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

### --- --- --- --- --- --- --- Project name --- --- --- --- --- --- --- ###

project(vector
        VERSION ${VERSION_NO}
        HOMEPAGE_URL "www.mozgan.me"
        DESCRIPTION "Library: vector"
        LANGUAGES "CXX"
)

enable_language(C)
enable_language(CXX)

### --- --- --- --- --- --- --- File globbing --- --- --- --- --- --- --- ###
file(GLOB_RECURSE headers include/*.h)
file(GLOB_RECURSE sources src/*.cpp)
set(SOURCES ${headers} ${sources})

### --- --- --- Include directory --- --- --- ###
include_directories(include)

### --- --- --- --- --- --- --- Compiler --- --- --- --- --- --- --- ###
### --- --- --- Compiler definition --- --- --- ###
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD 23)

### --- --- --- Compiler flags --- --- --- ###
### Compiler flags ###
# enable use of extra debugging information.
add_definitions(-g)

# enable all warnings about constructions that some users consider questionable.
add_definitions(-Wall)

# enable some extra warning flags that are not enabled by -Wall.
add_definitions(-Wextra)

# issue all the warnings demanded by strict ISO C/C++; reject all programs that
# use forbidden extensions, and some other programs that do not follow ISO C/C++.
add_definitions(-Wpedantic)

# do not warn about uses of functions, variables, and types marked as deprecated
# by using the deprecated attributes.
add_definitions(-Wno-deprecated-declarations)

# warn if declared variable, function, parameter etc. is not used.
add_definitions(-Wunused)

# warn when the order of member initializers given in the code does not match
# the order in which they must be executed.
add_definitions(-Wno-reorder)

# warn if the return type of a function has a type qualifier such as const.
add_definitions(-Wno-ignored-qualifiers)

# Warn about violations of the following style guidelines from Scott Meyers’
# Effective C++ series of books.
add_definitions(-Weffc++)

# Optimize even more.
#add_definitions(-O2)

# run the standard link-time optimizer. When invoked with source code,
# it generates GIMPLE and writes it to special ELF sections in the object file.
#add_definitions(-flto)

# During the link-time optimization, do not warn about type mismatches
# in global declarations from different compilation units.
#add_definitions(-Wlto-type-mismatch)

### --- --- --- --- --- --- Required libraries --- --- --- --- --- --- ###
### --- --- --- packages --- --- --- ###
find_package(Threads REQUIRED)
if (THREADS_FOUND)
    message(STATUS "CPP Threads found!")
else()
    message(FATAL_ERROR "CPP Threads not found!")
endif()

### --- --- --- static/shared libraries --- --- --- ###
#find_library(LIB_XXXXXX_FOUND
#    NAMES lib_xxxxxx       # lib_xxxxxx.a lib_xxxxxx.so lib_xxxxxx.dylib
#    HINTS /path/to/library
#)
#if (LIB_XXXXXX_FOUND)
#    message(STATUS "Library XXXXXX found!")
#else()
#    message(FATAL_ERROR "Library XXXXXX not found!")
#endif()

### --- --- --- --- --- --- Target definitions --- --- --- --- --- --- ###

### --- --- --- library --- --- --- ###
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
add_library(${CMAKE_PROJECT_NAME} SHARED ${SOURCES})

### --- --- --- --- --- --- Gtest --- --- --- --- --- --- ###
find_library(GTEST_LIBRARY gtest)
if (GTEST_LIBRARY)
    include_directories(test)
    enable_testing()
    add_subdirectory(test)
endif()

### --- --- --- --- --- --- Installation --- --- --- --- --- --- ###
#FILE(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
install(FILES ${headers} DESTINATION include COMPONENT headers)

set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE 1)
install(TARGETS ${PROJECT_NAME} DESTINATION lib COMPONENT libraries)

### --- --- --- --- --- --- Dependencies --- --- --- --- --- --- ###
# Warning: use the extension ".dylib" in macOS
#target_link_libraries(${CMAKE_PROJECT_NAME} lib_xxxxxx.dylib)
# Warning: in Linux no extension requested
#target_link_libraries(${CMAKE_PROJECT_NAME} add_some_lib)
#target_link_libraries(${CMAKE_PROJECT_NAME} pthread)

### --- --- --- --- --- --- Cpack --- --- --- --- --- --- ###
set(CPACK_PACKAGE_NAME lib${PROJECT_NAME}-dev)

set(CPACK_PACKAGE_VENDOR "Mehmet Ozgan")
set(CPACK_PACKAGE_CONTACT "mozgan@gmail.com")
set(CPACK_PACKAGE_HOMEPAGE_URL "www.mozgan.me")
set(CPACK_PACKAGE_VERSION ${VERSION_NO})
set(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${VERSION_NO}-${CMAKE_SYSTEM_PROCESSOR})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${CMAKE_PROJECT_DESCRIPTION}")

set(CPACK_COMPONENTS_ALL libraries headers)
set(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "${CMAKE_PROJECT_NAME}")
set(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "${CMAKE_PROJECT_DESCRIPTION}")
set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ Headers")
set(CPACK_COMPONENT_HEADERS_DESCRIPTION "C++ Headers for the library ${CMAKE_PROJECT_NAME}")

set(CPACK_PACKAGING_INSTALL_PREFIX "/usr/local")

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt DESTINATION /tmp)
endif()

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/README.txt DESTINATION /tmp)
endif()

set(CPACK_SET_DESTDIR ON)
set(CPACK_STRIP_FILES ON)
set(CPACK_SOURCE_STRIP_FILES ON)

if (${CMAKE_SYSTEM_NAME} MATCHES Darwin)
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/WELCOME.txt")
        set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_CURRENT_SOURCE_DIR}/WELCOME.txt")
    endif()

    set(CPACK_GENERATOR "productbuild;TGZ")
    set(CPACK_SOURCE_GENERATOR "TGZ")
elseif (${CMAKE_SYSTEM_NAME} MATCHES Linux)
    file(STRINGS "/etc/os-release" os_name REGEX "^NAME=")

    # create deb package for Debian and Ubuntu
    if(${os_name} MATCHES "^(NAME)=\"(Debian GNU/Linux|Ubuntu)\"$")
        set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${CMAKE_PROJECT_DESCRIPTION}")
        set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
        set(CPACK_DEBIAN_PACKAGE_DEPENDS "")
        set(CPACK_DEBIAN_PACKAGE_MAINTAINER ${CPACK_PACKAGE_VENDOR})
        set(CPACK_DEBIAN_PACKAGE_SECTION contrib/devel)
        set(CPACK_GENERATOR "DEB;TGZ")
    else()
        set(CPACK_GENERATOR "TGZ")
    endif()

    set(CPACK_SOURCE_GENERATOR "TGZ")
elseif (${CMAKE_SYSTEM_NAME} MATCHES FreeBSD)
    set(CPACK_FREEBSD_PACKAGE_DESCRIPTION "${CMAKE_PROJECT_DESCRIPTION}")
    set(CPACK_FREEBSD_PACKAGE_LICENSE "MIT")
    set(CPACK_FREEBSD_PACKAGE_DEPS "")
    set(CPACK_FREEBSD_PACKAGE_MAINTAINER ${CPACK_PACKAGE_VENDOR})
    set(CPACK_SOURCE_GENERATOR "TGZ")
else()
    set(CPACK_SOURCE_GENERATOR "TGZ")
endif()

set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/build/;${CMAKE_SOURCE_DIR}/_build/;${CMAKE_SOURCE_DIR}/.git/")

include(CPack)

