cmake_minimum_required(VERSION 3.16)
project(msp_calc)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(msp_calc main.cc)
target_link_libraries(msp_calc ${Geant4_LIBRARIES})

add_executable(range_calc range_calc.cc)
target_link_libraries(range_calc ${Geant4_LIBRARIES})

add_executable(sim_msp_range sim_msp_range.cc)
target_link_libraries(sim_msp_range ${Geant4_LIBRARIES})

# Copy file nang luong mac dinh vao thu muc build de chay ngay ma khong can
# tro duong dan tuyet doi
configure_file(${CMAKE_SOURCE_DIR}/energies_kev.txt
               ${CMAKE_BINARY_DIR}/energies_kev.txt COPYONLY)
