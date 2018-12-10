CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
 ADD_PARAVIEW_PLUGIN(DiffusionToConductivityTensor "1.0" 
   SERVER_MANAGER_XML vtkDiffusionToConductivityTensorFilter_SM.xml
   SERVER_MANAGER_SOURCES vtkDiffusionToConductivityTensorFilter.cxx
   REQUIRED_ON_SERVER)

# let cmake decide which c++-version is needed to satisfy all the features used in the source-code
target_compile_features(DiffusionToConductivityTensor PRIVATE cxx_range_for)
