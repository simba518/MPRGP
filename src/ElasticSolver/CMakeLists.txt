SET(LIBRARY_NAME elasticsolver)

AUX_SOURCE_DIRECTORY(. elasticsolver_src)

ADD_LIBRARY(${LIBRARY_NAME} ${elasticsolver_src})