project(tests LANGUAGES CXX)

aux_source_directory(. SOURCES_LIST)
foreach(filename ${SOURCES_LIST})
  get_filename_component(target ${filename} NAME_WE)
  add_executable(${target} ${filename})
  # link the different lib to exe
  target_link_libraries(${target} libsimu disSichel libcalc)
  target_compile_features(${target} PRIVATE cxx_std_17)
  # create test for each exe
  catch_discover_tests(${target} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/func_tests)
endforeach(filename)