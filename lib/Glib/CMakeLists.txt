project(libcalc)

aux_source_directory(. SOURCES_LIST)
# static lib are include in the final bin
add_library(${PROJECT_NAME} STATIC ${SOURCES_LIST})

SET(MY_DEBUG_OPTIONS "-Wall" "-Werror" "-Wfatal-errors" "-g")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:DEBUG>:${MY_DEBUG_OPTIONS}>")