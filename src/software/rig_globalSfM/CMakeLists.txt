#Add executable that computes:
# - features and descriptors
# - putatives matches
# - geometric filtered matches
#
ADD_EXECUTABLE(openMVG_main_computeMatchesRig main_computeMatchesRig.cpp)
TARGET_LINK_LIBRARIES(openMVG_main_computeMatchesRig
                     ${OpenMVG_LIBS})
SET_PROPERTY(TARGET openMVG_main_computeMatchesRig PROPERTY FOLDER OpenMVG/software)

INSTALL(TARGETS openMVG_main_computeMatchesRig DESTINATION bin/)

#Add library and executable that perform SfM:
# - the required input is matches computed by openMVG_main_computeMatches
#

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

INSTALL(TARGETS openMVG_SfM_Global_Rigid DESTINATION lib/)
INSTALL(TARGETS openMVG_main_rigidRig_GlobalSfM DESTINATION bin/)
