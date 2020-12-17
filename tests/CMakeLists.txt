project(tests LANGUAGES CXX)

#Toutes ces finctions sont des fonctions CMake sauf catch_discover_tests qui est implémenté dans catch.cmake dans lib/catch
aux_source_directory(. SOURCES_LIST) # prend le resultat de ls . (tous les fichiers du dossiers) et le met dans la variable SOURCE_LIST
foreach(filename ${SOURCES_LIST}) # boucle sur les noms de fichiers
  get_filename_component(target ${filename} NAME_WE) # enleve le .cpp du nom et  créé la variable target
  add_executable(${target} ${filename}) # créé ce qu'il faut pour faire un executable a partir de ce nom (créé un make avec la commande gcc...)
  # link the different lib to exe (comme les libs sont en static, il integre les .o/.a dans les exe)
  target_link_libraries(${target} libsimu disSichel libcalc) # créé ce qu'il faut pour linker
  target_compile_features(${target} PRIVATE cxx_std_17) # ajoute les options a mettre dans le make pour compiler
  # create test for each exe
  # fonction CMake faite par catch pour considérer chaque sous partie des tests (chaque TEST_ CASE) dans le framework de test (Ctest)
  catch_discover_tests(${target} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
endforeach(filename)
