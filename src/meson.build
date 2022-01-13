harmony_src_include_directories = include_directories('.')

harmony_config = configuration_data()
harmony_config.set('HARMONY_VERSION', meson.project_version())

harmony_gen_headers = [
  vcs_tag(
    input : 'LibraryGitHash.hpp.in',
    output : 'LibraryGitHash.hpp',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@HARMONY_GIT_SHA1@'),
  configure_file(
    input : 'LibraryVersion.hpp.in',
    output : 'LibraryVersion.hpp',
    configuration : harmony_config),
]

# sources + executable
harmony_main = executable(
  'harmony',
  files([
    'HarmonySettings.cpp',
    'LibraryInfo.cpp',
    'main.cpp',
    'SimpleBamParser.cpp',
  ]) + harmony_gen_headers,
  install : true,
  dependencies : harmony_lib_deps,
  include_directories : harmony_src_include_directories,
  cpp_args : harmony_flags)
