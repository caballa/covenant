set (USED_LIBS 
  avy
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(covenant covenant.cpp)
target_link_libraries (covenant ${USED_LIBS})
