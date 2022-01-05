add_library(xdr 
			STATIC)

target_sources(xdr
	PRIVATE
	xdrfile.c
	xdrfile_xtc.c
	xdrfile_trr.c
	xtc_seek.c
	trr_seek.c
	PUBLIC
	xdrfile.h
	xdrfile_xtc.h
	xdrfile_trr.h
	xtc_seek.h
	trr_seek.h
)

target_include_directories(xdr
	PUBLIC
	${CMAKE_CURRENT_LIST_DIR}
)

