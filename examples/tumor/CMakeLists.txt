add_example(
    tumor
    GALOIS
    SRC
        vasculature/plot.cpp
        vasculature/vasculature.cpp
        tumor.cpp
        tumor_main.cpp
)
add_example(
    tumor_3d
    GALOIS
    SRC #
        3d/vasculature_parser.cpp
        3d/main.cpp
)
