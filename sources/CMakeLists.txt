project(libsimu)

aux_source_directory(. SOURCES_LIST)
# static lib are include in the final bin
add_library(${PROJECT_NAME} STATIC ${SOURCES_LIST})

SET(MY_DEBUG_OPTIONS "-Wall" "-Werror" "-Wfatal-errors" "-ansi" "-pedantic")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:DEBUG>:${MY_DEBUG_OPTIONS}>")

SET(MY_RELEASE_WITH_DEBUG_OPTIONS "-Wall" "-Werror" "-Wfatal-errors")
target_compile_options(${PROJECT_NAME} PUBLIC "$<$<CONFIG:RELWITHDEBINFO>:${MY_RELEASE_WITH_DEBUG_OPTIONS}>")