add_custom_target(CopyData        
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/data
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data        
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )

add_dependencies(LoveNumbers CopyData)        

add_executable(ex1 ex1.cpp)
target_link_libraries(ex1 LoveNumbers)

add_executable(ex2 ex2.cpp)
target_link_libraries(ex2 LoveNumbers)