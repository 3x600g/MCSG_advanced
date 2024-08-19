cmake_minimum_required(VERSION 3.23)
project(MCSG1_2)

set(CMAKE_CXX_STANDARD 14)

add_executable(MCSG1_2
    main.cpp
    Surface.cpp
    Material.cpp
    CellCard.cpp
    Calculation.cpp
    Utils.cpp
)
