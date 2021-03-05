project(tskit)

aux_source_directory(. SOURCES_LIST)
# static lib are include in the final bin
add_library(${PROJECT_NAME} STATIC ${SOURCES_LIST})