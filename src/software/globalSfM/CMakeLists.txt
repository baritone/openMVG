
ADD_LIBRARY(openMVG_SfM_Global SfMGlobalEngine.cpp)

ADD_EXECUTABLE(openMVG_main_GlobalSfM main_GlobalSfM.cpp)
TARGET_LINK_LIBRARIES(
  openMVG_main_GlobalSfM
  openMVG_SfM_Global
  ${OpenMVG_LIBS})

SET_PROPERTY(TARGET openMVG_main_GlobalSfM PROPERTY FOLDER OpenMVG/software)
SET_PROPERTY(TARGET openMVG_SfM_Global PROPERTY FOLDER OpenMVG/OpenMVG)

IF(APPLE)
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/third_party/opengv/include/ ${CMAKE_SOURCE_DIR}/third_party/opengv/third_party/)
ELSE()
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/third_party/opengv/include/)
ENDIF()

FIND_LIBRARY(OPENGV_LIBRARY opengv ${CMAKE_BINARY_DIR}/third_party/opengv/lib)
FIND_LIBRARY(OPENGV_LIBRARY random_generators ${CMAKE_BINARY_DIR}/third_party/opengv/lib)

ADD_LIBRARY(openMVG_SfM_Global_Rigid SfMRigidGlobalEngine.cpp)

ADD_EXECUTABLE(openMVG_main_rigidRig_GlobalSfM main_rigidRig_GlobalSfM.cpp)
TARGET_LINK_LIBRARIES(
  openMVG_main_rigidRig_GlobalSfM
  openMVG_SfM_Global_Rigid
  opengv
  random_generators
  ${OpenMVG_LIBS})

SET_PROPERTY(TARGET openMVG_main_rigidRig_GlobalSfM PROPERTY FOLDER OpenMVG/software)
SET_PROPERTY(TARGET openMVG_SfM_Global_Rigid PROPERTY FOLDER OpenMVG/OpenMVG)
