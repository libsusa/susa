add_test (NAME types COMMAND types)
add_executable (types types.cpp)
target_link_libraries (types susa)

add_test (NAME base COMMAND base)
add_executable (base base.cpp)
target_link_libraries (base susa)

add_test (NAME linalg COMMAND linalg)
add_executable (linalg linalg.cpp)
target_link_libraries (linalg susa)
