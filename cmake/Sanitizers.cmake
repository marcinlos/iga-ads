# Enable address sanitizer
target_compile_options(project-options
  INTERFACE
  $<$<CONFIG:Debug>:-fsanitize=address,undefined;-fno-omit-frame-pointer>)

target_link_options(project-options
  INTERFACE
  $<$<CONFIG:Debug>:-fsanitize=address,undefined>)

