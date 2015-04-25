find_program(PYTHON "python")

option (USE_PY_SETUP "Use python distutils to install python modules" OFF)

if (PYTHON AND NOT USE_PY_SETUP)
  install (PROGRAMS brunch.py    DESTINATION bin )
  install (PROGRAMS dotToPs.py DESTINATION bin )
  install (PROGRAMS covenant-par DESTINATION bin )
  install (FILES LogManager.py   DESTINATION bin)
endif()

if (PYTHON AND USE_PY_SETUP)
  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

  set(DEPS        
    brunch.py
    dotToPs.py
    covenant-par
    LogManager.py
    ${SETUP_PY_IN})
  
  configure_file(${SETUP_PY_IN} ${SETUP_PY})

  add_custom_command(OUTPUT ${OUTPUT}
    COMMAND ${PYTHON} ${SETUP_PY} build
    COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
    DEPENDS ${DEPS})

  add_custom_target(python ALL DEPENDS ${OUTPUT})

  install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install --install-lib ${CMAKE_INSTALL_PREFIX}/bin --install-scripts ${CMAKE_INSTALL_PREFIX}/bin)")
endif()
