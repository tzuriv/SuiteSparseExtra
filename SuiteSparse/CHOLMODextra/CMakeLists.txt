PROJECT(CHOLMODextra)

include_directories("${suitesparseconfig_SOURCE_DIR}")
include_directories("${AMD_SOURCE_DIR}/Include")
include_directories("${CAMD_SOURCE_DIR}/Include")
include_directories("${COLAMD_SOURCE_DIR}/Include")
include_directories("${CCOLAMD_SOURCE_DIR}/Include")
include_directories("${CHOLMOD_SOURCE_DIR}/Include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# Headers:
include_directories("${CHOLMODextra_SOURCE_DIR}/Include")
FILE(GLOB LIBHDRS "${CHOLMODextra_SOURCE_DIR}/Include/*.h")

FILE(GLOB LIBSRCS "Source/*.c")

declare_suitesparse_library(cholmod-extra "${LIBSRCS}" "${LIBHDRS}"
	TARGET_PRIVATE_LINK	${SuiteSparse_LINKER_LAPACK_BLAS_LIBS}
	TARGET_PUBLIC_LINK 	${SuiteSparse_LINKER_METIS_LIBS} amd camd colamd ccolamd
)

