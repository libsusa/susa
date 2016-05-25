add_executable (simple_matrix simple_matrix.cpp)
target_link_libraries (simple_matrix susa)

add_executable (bcjr_bpsk_cc75 bcjr_bpsk_cc75.cpp)
target_link_libraries (bcjr_bpsk_cc75 susa)

add_executable (qam_ser_gaussian qam_ser_gaussian.cpp)
target_link_libraries (qam_ser_gaussian susa)

add_executable (dijkstra dijkstra.cpp)
target_link_libraries (dijkstra susa)

