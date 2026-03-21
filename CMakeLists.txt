cmake_minimum_required(VERSION 3.10)
project(LinearProgramming)

# Устанавливаем стандарт C++11 (используется в коде)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Добавляем флаги предупреждений для GCC/Clang (необязательно)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# Создаём исполняемый файл из всех .cpp файлов
add_executable(lp_solver
    main.cpp
    lp_utils.cpp
    simplex.cpp
    lp_types.h        # заголовки перечислены для удобства IDE
    lp_utils.h
    simplex.h
)

# Указываем директорию, где лежат заголовочные файлы (текущая)
target_include_directories(lp_solver PRIVATE .)
