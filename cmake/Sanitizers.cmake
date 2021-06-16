# Enable address sanitizer
target_compile_options(project_options
  INTERFACE
  $<$<CONFIG:Debug>:-fsanitize=address,undefined;-fno-omit-frame-pointer>)

target_link_options(project_options
  INTERFACE
  $<$<CONFIG:Debug>:-fsanitize=address,undefined>)

