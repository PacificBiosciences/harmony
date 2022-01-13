project(
  'harmony',
  'cpp',
  version : '0.1.0-99',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++20',
    'b_ndebug=if-release'],
  license : 'BSD-3',
  meson_version : '>= 0.57.0')

# CXXFLAGS
harmony_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Werror=strict-aliasing',
  '-Wno-unused-parameter',
  '-Wno-unused-local-typedef']
    if cpp.has_argument(cflag)
      harmony_flags += cflag
    endif
endforeach

# dependencies #
## threads
harmony_thread_dep = dependency('threads', required : true)
## boost
harmony_boost_dep = dependency('boost', required : true)
## zlib
harmony_zlib_dep = dependency('zlib', required : true, fallback : ['zlib', 'zlib_dep'])
## pbbam
harmony_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])
## pbcopper
harmony_pbcopper_dep = dependency('pbcopper', fallback : ['pbcopper', 'pbcopper_dep'])
# htslib
harmony_htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])

harmony_lib_deps = [
  harmony_thread_dep,
  harmony_boost_dep,
  harmony_zlib_dep,
  harmony_pbbam_dep,
  harmony_pbcopper_dep,
  harmony_htslib_dep]

subdir('src')
