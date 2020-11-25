add_test (NAME types COMMAND types)
add_executable (types types.cpp)
target_link_libraries (types susa)

add_test (NAME base COMMAND base)
add_executable (base base.cpp)
target_link_libraries (base susa)

add_test (NAME linalg COMMAND linalg)
add_executable (linalg linalg.cpp)
target_link_libraries (linalg susa)

add_test (NAME signal COMMAND signal)
add_executable (signal signal.cpp)
target_link_libraries (signal susa)

add_test (NAME channel COMMAND channel)
add_executable (channel channel.cpp)
target_link_libraries (channel susa)

add_test (NAME ccode COMMAND ccode)
add_executable (ccode ccode.cpp)
target_link_libraries (ccode susa)

add_test (NAME memory COMMAND memory)
add_executable (memory memory.cpp)
target_link_libraries (memory susa)

add_test (NAME fft COMMAND fft)
add_executable (fft fft.cpp)
target_link_libraries (fft susa)

add_test (NAME matrix COMMAND matrix)
add_executable (matrix matrix.cpp)
target_link_libraries (matrix susa)