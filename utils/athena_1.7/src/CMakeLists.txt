include_directories(${ATHENA_SOURCE_DIR}/src)


# Targets

# athena
add_executable(athena main.c 2pcf.c config.c errorlist.c gal_cat.c io.c kdtree.c mathstuff.c nrcomplex.c)
target_link_libraries(athena ${LIBS})

# venice
add_executable(venice venice.c)
target_link_libraries(venice gsl gslcblas ${LIBS})



# Installation
install(TARGETS athena venice DESTINATION bin)



