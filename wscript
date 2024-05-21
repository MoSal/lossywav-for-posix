#! /usr/bin/env python

# waf imports
from waflib.Configure import conf

def options(opt):
    opt.load('compiler_cxx')

    ins_gr = opt.get_option_group('Installation and uninstallation options')
    bld_gr = opt.get_option_group('Build and installation options')
    conf_gr = opt.get_option_group('Configuration options')
 
    def_enable_compiler_warnings = False
    conf_gr.add_option(
            '--enable-compiler-warnings',
            dest = 'ENABLE_COMPILER_WARNINGS',
            default = def_enable_compiler_warnings,
            action= "store_true",
            help = "Enable compiler warnings. (default: %s)" % def_enable_compiler_warnings
            )

    def_werror = False
    conf_gr.add_option(
            '--werror',
            dest = 'WERROR',
            default = def_werror,
            action= "store_true",
            help = "Consider warnings fatal. (default: %s)" % def_werror
            )

    def_disable_compile_opt = False
    conf_gr.add_option(
            '--disable-compile-optimizations',
            dest = 'DISABLE_CXXFLAGS_OPT',
            default = def_disable_compile_opt,
            action= "store_true",
            help = "Don't check/set compile optimization flags. (default: %s)" % def_disable_compile_opt
            )

    def_disable_link_opt = False
    conf_gr.add_option(
            '--disable-link-optimizations',
            dest = 'DISABLE_LINK_OPT',
            default = def_disable_link_opt,
            action= "store_true",
            help = "Don't check/set link optimization flags. (default: %s)" % def_disable_link_opt
            )

    def_disable_lto = False
    conf_gr.add_option(
            '--disable-lto',
            dest = 'DISABLE_LTO',
            default = def_disable_lto,
            action= "store_true",
            help = "Don't check/set lto flags. (default: %s)" % def_disable_lto
            )

    def_enable_debug = False
    conf_gr.add_option(
            '--enable-debug',
            dest = 'ENABLE_DEBUG',
            default = def_enable_debug,
            action= "store_true",
            help = "Set debug flags. (default: %s)" % def_enable_debug
            )

    def_enable_fftw3 = False
    conf_gr.add_option(
            '--enable-fftw3',
            dest = 'ENABLE_FFTW3',
            default = def_enable_fftw3,
            action= "store_true",
            help = "Compile and link against libfftw3. (default: %s)" % def_enable_fftw3
            )

    def_fftw3_cxxflags = None # Use pkg-config or env
    conf_gr.add_option(
            '--fftw3-cxxflags',
            dest = 'FFTW3_CXXFLAGS',
            default = def_fftw3_cxxflags,
            action= "store",
            help = "Skip pkg-config and set fftw3 cxxflags explicitly (default: %s)" % def_fftw3_cxxflags
            )

    def_fftw3_libs = None # Use pkg-config or env
    conf_gr.add_option(
            '--fftw3-libs',
            dest = 'FFTW3_LIBS',
            default = def_fftw3_libs,
            action= "store",
            help = "Skip pkg-config and set fftw3 libs explicitly (default: %s)" % def_fftw3_libs
            )

#------------------------------------------------------------------------------

@conf
def check_func(conf, f_name, h_name, mandatory=True):
    fragment = ''
    if len(h_name) > 0:
        fragment += '#include <' + h_name + '>\n'
    fragment +='int main() { void(*p)(void) = (void(*)(void)) ' + f_name + '; return !p; }\n'

    msg = 'Checking for ' + f_name + '()'
    define_name = 'HAVE_' + f_name.upper().replace('::', '_')

    conf.check_cxx(fragment=fragment, define_name=define_name,  mandatory=mandatory, msg=msg)

@conf
def check_api(conf):
    if conf.env['DEST_OS'] != 'win32':
        print('Checking API support:')
        check_func(conf, "std::chrono::steady_clock::now", "chrono")
        check_func(conf, "setpriority", "sys/resource.h")
        check_func(conf, "stat", "sys/stat.h")
        check_func(conf, "chmod", "sys/stat.h")
        check_func(conf, "nanosleep", "ctime")
        check_func(conf, "sincos", "math.h", False)

@conf
def check_flags(conf):
    # Load this before checking flags
    conf.load('compiler_cxx')

    # Set warning flags 1st so -Werror catches all warnings
    if conf.options.ENABLE_COMPILER_WARNINGS:
        check_warning_cxxflags(conf)

    check_required_flags(conf)

    if conf.options.ENABLE_DEBUG:
        check_debug_cxxflags(conf)
    else:
        if not conf.options.DISABLE_CXXFLAGS_OPT:
            check_opt_cxxflags(conf)
        if not conf.options.DISABLE_LTO:
            check_lto_flags(conf)

    if not conf.options.DISABLE_LINK_OPT:
        check_link_flags(conf)

@conf
def check_required_flags(conf):
    conf.check_cxx(cxxflags = '-std=c++11', uselib_store='LOSSYWAV_REQUIRED', mandatory=False)
    conf.check_cxx(cxxflags = '-O2', uselib_store='LOSSYWAV_REQUIRED', mandatory=False)
    conf.check_cxx(cxxflags = '-pipe', uselib_store='LOSSYWAV_REQUIRED', mandatory=False)

    conf.check_cxx(cxxflags = '-fPIE', uselib_store='LOSSYWAV_REQUIRED', mandatory=False)
    conf.check_cxx(linkflags = '-pie', uselib_store='LOSSYWAV_REQUIRED', mandatory=False)

    if conf.env['CXXFLAGS_LOSSYWAV_REQUIRED']:
        conf.env.append_value('CXXFLAGS', conf.env['CXXFLAGS_LOSSYWAV_REQUIRED'])

@conf
def check_warning_cxxflags(conf):
    print('Checking for warning CXXFLAGS support:')

    warn_flags = [
            ['-Wall'],
            ['-Wextra'],
            ['-Wmissing-format-attribute'],
    ]

    for w in warn_flags:
        conf.check_cxx(cxxflags = w, uselib_store='LOSSYWAV_WARNING', mandatory=False)

    # Set -Werror if requested
    if conf.options.WERROR:
        conf.check_cxx(cxxflags = ['-Werror'], uselib_store='LOSSYWAV_WARNING', mandatory=False)
    else:
        conf.check_cxx(cxxflags = ['-Wno-error'], uselib_store='LOSSYWAV_WARNING', mandatory=False)


    if conf.env['CXXFLAGS_LOSSYWAV_WARNING']:
        conf.env.append_value('CXXFLAGS', conf.env['CXXFLAGS_LOSSYWAV_WARNING'])
    else:
        conf.fatal('None of the warning CXXFLAGS are supported by the compiler!')

@conf
def check_debug_cxxflags(conf):
    print('Checking for debug CXXFLAGS support:')

    debug_flags = [
            ['-Og'],
            ['-ggdb'],
            ['-fvar-tracking-assignments'],
            ['-fno-omit-frame-pointer'],
            ['-fstack-protector-strong'],
            ['-fstack-check']
    ]

    for d in debug_flags:
        conf.check_cxx(cxxflags = d, uselib_store='LOSSYWAV_DEBUG', mandatory=False)

    if conf.env['CXXFLAGS_LOSSYWAV_DEBUG']:
        conf.env.append_value('CXXFLAGS', conf.env['CXXFLAGS_LOSSYWAV_DEBUG'])
    else:
        conf.fatal('None of the debug CXXFLAGS are supported by the compiler!')

@conf
def check_opt_cxxflags(conf):
    print('Checking for optimized CXXFLAGS support:')

    conf.check_cxx(cxxflags = '-Ofast', uselib_store='LOSSYWAV_OPT', mandatory=False)
    if not conf.env['CXXFLAGS_LOSSYWAV_OPT']:
        conf.check_cxx(cxxflags = '-O3', uselib_store='LOSSYWAV_OPT', mandatory=False)
        conf.check_cxx(cxxflags = '-ffast-math', uselib_store='LOSSYWAV_OPT', mandatory=False)

    if conf.env['CXXFLAGS_LOSSYWAV_OPT']:
        if '-O2' in conf.env['CXXFLAGS']:
            conf.env['CXXFLAGS'].remove('-O2')
        conf.env.append_value('CXXFLAGS', conf.env['CXXFLAGS_LOSSYWAV_OPT'])

@conf
def check_link_flags(conf):
    print('Checking for optimized LINKFLAGS support:')

    linkflags = [
            ['-Wl,-O1'],
            ['-Wl,--sort-common'],
            ['-Wl,--as-needed'],
            ['-Wl,-z,relro'],
            ['-Wl,-z,now'],
    ]

    for l in linkflags:
        conf.check_cxx(linkflags = l, uselib_store='LOSSYWAV', mandatory=False)

    if conf.env['LINKFLAGS_LOSSYWAV']:
        conf.env.append_value('LINKFLAGS', conf.env['LINKFLAGS_LOSSYWAV'])

@conf
def check_lto_flags(conf):
    print('Checking for LTO CXXFLAGS/LINKFLAGS support:')
    conf.check_cxx(cxxflags = '-flto', uselib_store='LOSSYWAV_LTO', mandatory=False)
    conf.check_cxx(linkflags = '-flto', uselib_store='LOSSYWAV_LTO', mandatory=False)

    if conf.env['CXXFLAGS_LOSSYWAV_LTO'] or conf.env['LINKFLAGS_LOSSYWAV_LTO']:
        if conf.env['CXXFLAGS_LOSSYWAV_LTO']:
            conf.env.append_value('CXXFLAGS', conf.env['CXXFLAGS_LOSSYWAV_LTO'])

        if conf.env['LINKFLAGS_LOSSYWAV_LTO']:
            conf.env.append_value('LINKFLAGS', conf.env['LINKFLAGS_LOSSYWAV_LTO'])
    else:
        print('lto flags not supported by the compiler!')


@conf
def check_pkg_deps(conf):
    print('Check dependencies:')

    # Make sure these exist exist
    for v in ['INCLUDES', 'RPATH', 'CXXFLAGS', 'LDFLAGS' 'LIB' 'LIBPATH']:
        if not v in conf.env:
            conf.env[v] = []

    if conf.options.ENABLE_FFTW3:
        check_fftw3(conf)

@conf
def check_fftw3(conf):
    pkg_name = 'fftw3'
    check_args = ['--cflags', '--libs']
    min_ver = '3.3.0'
    check_pkg(conf, pkg_name, check_args, min_ver)

@conf
def check_pkg(conf, pkg_name, check_args, min_ver):

    conf_opts_dict = eval( str(conf.options) )

    opt_cxxflags_var = pkg_name.upper() + '_CXXFLAGS'
    opt_libs_var = pkg_name.upper() + '_LIBS'

    opt_cxxflags = (opt_cxxflags_var in conf_opts_dict) and  conf_opts_dict[opt_cxxflags_var] != None
    opt_libs = (opt_libs_var in conf_opts_dict) and conf_opts_dict[opt_libs_var] != None

    conf.start_msg('Checking %s:' % pkg_name)

    if opt_cxxflags and opt_libs:
        conf.end_msg('user-provided')
        conf.env['CXXFLAGS'] += conf_opts_dict[opt_cxxflags_var].split(' ')
        conf.env['LDFLAGS'] += conf_opts_dict[opt_libs_var].split(' ')
    elif opt_cxxflags or opt_libs:
        conf.fatal('Either set both %s and %s or let pkg-config do the checking.' % (opt_cxxflags_var, opt_libs_var))
    else:
        conf.end_msg('pkg-config')
        conf.check_cfg(package = pkg_name, variables = ['includedir', 'prefix'])
        # Store without version
        conf.check_cfg(package = pkg_name + ' >= ' + min_ver, args = check_args, uselib_store=pkg_name.upper())

        defines_var = 'DEFINES_' + pkg_name.upper()
        if conf.env[defines_var]:
            conf.env.DEFINES += conf.env[defines_var]

        includes_var = 'INCLUDES_' + pkg_name.upper()
        if conf.env[includes_var]:
            conf.env.INCLUDES += conf.env[includes_var]

        if conf.env[pkg_name.upper() + '_includedir']:
            conf.env.INCLUDES += [ conf.env[pkg_name.upper() + '_includedir'] ]

        # NetBSD relies on RPATH instead of ldconfig
        # This would only be set if -Wl,-R<dir> (or -Wl,-rpath<DIR>) is a part of Libs
        rpath_var = 'RPATH_' + pkg_name.upper()
        if conf.env[rpath_var]:
            conf.env.RPATH += conf.env[rpath_var]

        libpath_var = 'LIBPATH_' + pkg_name.upper()
        if conf.env[libpath_var]:
            conf.env.LIBPATH += conf.env[libpath_var]

        cxxflags_var = 'CXXFLAGS_' + pkg_name.upper()
        if conf.env[cxxflags_var]:
            conf.env.CXXFLAGS += conf.env[cxxflags_var]

        lib_var = 'LIB_' + pkg_name.upper()
        if conf.env[lib_var]:
            conf.env.LIB += conf.env[lib_var]

#------------------------------------------------------------------------------

def configure(conf):
    check_flags(conf)
    check_api(conf)
    check_pkg_deps(conf)

#------------------------------------------------------------------------------

def build(bld):

    bld.objects(
            source = [
                'units/fftw_interface.cpp',
                'units/nCore.cpp',
                'units/nFFT.cpp',
                'units/nFillFFT.cpp',
                'units/nInitialise.cpp',
                'units/nOutput.cpp',
                'units/nParameter.cpp',
                'units/nProcess.cpp',
                'units/nRemoveBits.cpp',
                'units/nSGNS.cpp',
                'units/nShiftBlocks.cpp',
                'units/nSpreading.cpp',
                'units/nWav.cpp',
                ],
            target = ['lossywav-objs']
            )

    bld.program(
            use = ['lossywav-objs'],
            source = ['lossyWAV.cpp'],
            target = 'lossywav'
            )

#------------------------------------------------------------------------------
