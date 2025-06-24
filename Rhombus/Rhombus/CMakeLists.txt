set(RHOMBUS_SRC_LIST
    "status.cpp"
    "seal_api.cpp"
    "matvec.cpp"
    "matrix.cpp"
    "matmul.cpp"
    "matmul_protocol.cpp"
)

add_library(Rhombus ${RHOMBUS_SRC_LIST})
target_include_directories(Rhombus PUBLIC ${SCI_INCLUDE_DIRS})
target_link_libraries(Rhombus PUBLIC SCI-common SEAL::seal)
target_compile_definitions(Rhombus PUBLIC SCI_OT=1 USE_CHEETAH=1)